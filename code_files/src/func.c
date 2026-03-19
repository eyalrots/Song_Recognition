#include "../include/func.h"
#include <stdint.h>
#include <stdio.h>

int read_file(const char *input_file, unsigned char** out, long* size, wav_header_t* header) {
    unsigned char *data_buffer = NULL;
    FILE* wav_file = fopen(input_file, "rb");

    if (wav_file==NULL) {
        perror("Error opening file\n");
        return 1;
    }

    size_t bytes_read = fread(header, 1, sizeof(wav_header_t), wav_file);
    if (bytes_read < sizeof(wav_header_t)) {
        fprintf(stderr, "Error: Could not read WAV header.\n");
        fclose(wav_file);
        return 1;
    }

    // Validate WAV header
    if (strncmp(header->riff_header, "RIFF", 4) != 0 || 
        strncmp(header->wave_header, "WAVE", 4) != 0) {
            fprintf(stderr, "Error: Not a valid RIFF/WAVE file.\n");
            fclose(wav_file);
            return 1;
    }

    char chunk_header[4];
    uint32_t chunk_size;

    while (1) {
        // Reading chunk header and size
         if (fread(chunk_header, 1, 4, wav_file) != 4) {
            fprintf(stderr, "Error: Reached end of file without finding 'data' chunk.\n");
            fclose(wav_file);
            return 1;
        }
        if (fread(&chunk_size, 1, 4, wav_file) != 4) {
            fprintf(stderr, "Error: Could not read chunk size.\n");
            fclose(wav_file);
            return 1;
        }

        // Check if this is the 'data' chunk
        if (strncmp(chunk_header, "data", 4) == 0) {
            printf("Found 'data' chunk. Size: %u bytes.\n", chunk_size);
            break; // Found it, exit the loop
        }

        // If not 'data', skip this chunk's data to get to the next chunk header
        printf("Skipping '%.4s' chunk of size %u.\n", chunk_header, chunk_size);
        fseek(wav_file, chunk_size, SEEK_CUR);
    }

    // Print WAV file information
    printf("--- WAV File Header ---\n");
    printf("Sample Rate:        %u Hz\n", header->sample_rate);
    printf("Bits Per Sample:    %u bits\n", header->bits_per_sample);
    printf("Channels:           %u\n", header->num_channels);
    printf("Audio Format:       %u (1=PCM)\n", header->audio_format);
    printf("-----------------------\n\n");

    // Allocate a buffer to hold the audio data
    data_buffer = calloc(chunk_size, sizeof(char));
    if (!data_buffer) {
        perror("Error allocating memory for data\n");
        fclose(wav_file);
        return 1;
    }

    // Read audio data
    bytes_read = fread(data_buffer, 1, chunk_size, wav_file);
    if (bytes_read < chunk_size) {
        fprintf(stderr, "Error reading audio data.\n");
    } else {
        printf("Successfully read %zu bytes of audio data.\n", bytes_read);   
    }

    *out = data_buffer;
    *size = chunk_size;

    return 0;
}

double get_sample(unsigned char* samples, int bytes_pet_sample) {
    int i=0;
    uint32_t sample = 0;
    for (i=bytes_pet_sample-1; i>=0; i--) {
        sample = (sample << 8) | samples[i];
    }
    return (double)sample;
}

int convert_to_samples(unsigned char* data_buffer, double** sample_data, long* size, int bits_per_sample) {
    int i=0;
    int bytes_per_sample = bits_per_sample / 8;
    int16_t raw = 0;
    int byte_idx = 0;

    if (!data_buffer) {
        return 1;
    }

    *size = *size / bytes_per_sample;
    *sample_data = (double*)calloc(*size, sizeof(double));
    if (!(*sample_data)) return -1;
    for (i=0; i < *size - 1; i++) {
        byte_idx = i * bytes_per_sample;
        raw = (int16_t)(data_buffer[byte_idx] | (data_buffer[byte_idx+1] << 8));
        //(*sample_data)[i] = get_sample(&data_buffer[i*bytes_per_sample], bytes_per_sample);
        (*sample_data)[i] = (double)raw;
    }

    return 0;
}

int hamming(double* sub_window, int sub_size) {
    double hamming_value = 0;
    int i = 0;

    if (!sub_window || !sub_size) {
        return 1;
    }

    for (i=0; i < sub_size; i++) {
        hamming_value = 0.54 - 0.46 * cos(2 * M_PI * i / (sub_size-1));
        sub_window[i] *= hamming_value;
    }

    return 0;
}

int create_periodogram(uint32_t fft_size, const fftw_complex* out, double* periodogram) {
    /* Receives a output of an FFT.
     * Calculates the squared magnitude of said FFT.
     */
    int i = 0;
    double real = 0;
    double imag = 0;
    double mag_sq = 0;

    if (!out || !periodogram) {
        return -1;
    }
    /* Calculate */
    for (i = 0; i < (int)fft_size; i++) {
        real = out[i][0];
        imag = out[i][1];

        mag_sq = (real * real) + (imag * imag);

        /* Multiply by two for non DC frequencies */
        if (i) mag_sq *= 2.0;

        periodogram[i] = mag_sq;
    }
    return 0;
}

int fft_to_mel(const uint32_t fft_size, double* mel_bins, int num_of_mel_bins) {
    double min_mel = 0;
    double max_mel = 0;
    double jump = 0;
    int i = 0;

    if (!mel_bins || !num_of_mel_bins) {
        return -1;
    }
    // map Hz to Mel
    min_mel = 2595 * log(1 + ((double)MIN_FREQ / 700));
    max_mel = 2595 * log(1 + ((double)MAX_FREQ / 700));

    // create points on mel scale spaced evenly
    jump = (max_mel - min_mel) / num_of_mel_bins;
    mel_bins[0] = min_mel;
    for (i = 1; i < num_of_mel_bins; i++) {
        mel_bins[i] = mel_bins[i-1] + jump;
    }

    // convert back to Hz
    for (i = 0; i < num_of_mel_bins+1; i++) {
        mel_bins[i] = (int)(700 * exp(((double)mel_bins[i] / 2595) - 1));
    }

    // Convert to idx in original array -> idx = f * size / f_s(=sample rate)
    for (i = 0; i < num_of_mel_bins; i++) {
        mel_bins[i] = mel_bins[i] * fft_size / MAX_FREQ;
    }

    return 0;
}

int threshold_function(double* psd, double* a_out, double* b_out) {
/**
 * Calculates a(t_n) and b(t_n) based on the provided formulas.
 * B: Number of bins (summation limit)
 * periodogram_power: Array of power values P(t_n)[i]
 * a_out, b_out: Pointers to store the results
 */
    double sum_lnP = 0;
    double sum_j = 0;
    double sum_j_sq = 0;
    double sum_j_lnP = 0;
    double lnP = 0;
    double denominator;
    int N = LINEKEY_SIZE - 1;

    if (!psd || !a_out || !b_out) {
        fprintf(stderr, "Error: could not resolve pointers for inputs to threshold-function.\n");
        return -1;
    }

    for (int i = 2; i <= LINEKEY_SIZE; i++) {
        int j = i - 1; // This represents (i-1) in the formula
        
        /* P(t_n)[i] - 1-based indexing in formula, 
         * so we access index i-1 in the C array.
         */
        
        // Use a small epsilon to avoid log(0)
        lnP = log(psd[i-1] + 1e-12);

        sum_lnP += lnP;
        sum_j += (double)j;
        sum_j_sq += (double)(j * j);
        sum_j_lnP += (double)j * lnP;
    }

    // Common denominator for both formulas
    denominator = (N * sum_j_sq) - (sum_j * sum_j);

    if (fabs(denominator) < 1e-15) {
        *a_out = 0.0;
        *b_out = 0.0;
        return 0;
    }

    // Formula for a(t_n)
    *a_out = (sum_lnP * sum_j_sq - sum_j * sum_j_lnP) / denominator;

    /* Formula for b(t_n) 
     * Note: In the paper there is an 'n' in the denominator for b(t_n), 
     * but typically in linear regression this is B. 
     * I will use B as per standard least-squares derivation.
     */
    *b_out = (N * sum_j_lnP - sum_j * sum_lnP) / denominator;

    return 0;
}

int new_linekey_generation(double* window_data, int sigma_size, int sub_size, int sample_rate, uint64_t* out_linekey) {
    (void)sample_rate; /* Unused. */
    /* General variables */
    int return_val = 0;
    int i = 0;
    int m = 0;
    int k = 0;
    /* FFT variables */
    double in[sub_size];
    fftw_complex *out = NULL;
    fftw_plan fft_plan;
    uint32_t fft_size = sub_size / 2 + 1;

    /* Periodograms variables */
    uint32_t step_size = 0;
    double** periodograms = NULL; // An array of periodograms
    double* window_ptr = NULL;
    
    /* Mel variables */
    double after_mel[NUM_OF_PERIODOGRAMS][LINEKEY_SIZE];
    double mel_bins[LINEKEY_SIZE];

    /* PSD variables */
    double psd[LINEKEY_SIZE];

    /* Threshold variables */
    double a = 0;
    double b = 0;
    double cur_threshold = 0;
    int linekey_bits[LINEKEY_SIZE];

    if (!window_data || !sigma_size || !sub_size || !out_linekey) {
        fprintf(stderr, "Error: Input values are invalid. Could not generate LineKey.\n");
        return_val = -1;
        goto out;
    }

    /* Prepare FFT structures */
    out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * fft_size);
    if (!out) {
        fprintf(stderr, "Error: Failed to allocate memory for FFT.\n");
        return_val = -1;
        goto out;
    }
    /* Each time we need to calculate:
     * 1. Put the input samples in -in- array.
     * 2. fftw_execute(fft_plan);
     */
    fft_plan = fftw_plan_dft_r2c_1d(sub_size, in, out, FFTW_ESTIMATE);

    /* Create Periodograms */
    periodograms = (double **) calloc(NUM_OF_PERIODOGRAMS, sizeof(double *));
    if (!periodograms) {
        fprintf(stderr, "Error: Could not allocate memory for periodograms.\n");
        return_val = -1;
        goto out;
    }
    memset(in, 0, sizeof(in));

    step_size = (sigma_size - sub_size) / NUM_OF_PERIODOGRAMS;
    window_ptr = window_data;
    
    for (i = 0; i < NUM_OF_PERIODOGRAMS; i++) {
        /* Set current window */
        if (i) window_ptr += step_size;
        memcpy(in, window_ptr, sizeof(in));
        /* Allocate space for new periodogram */
        periodograms[i] = (double *) calloc(fft_size, sizeof(double));
        /* Calculations */
        hamming(in, sub_size);
        fftw_execute(fft_plan);
        if (create_periodogram(fft_size, out, periodograms[i]) < 0) {
            return_val = -1;
            fprintf(stderr, "Error: Failed to create priodograms.\n");
            goto out;
        }
    }
    
    /* Mel Filters */
    memset(after_mel, 0, sizeof(after_mel));
    memset(mel_bins, 0, sizeof(mel_bins));
    fft_to_mel(fft_size, mel_bins, LINEKEY_SIZE);

    for (k = 0; k < NUM_OF_PERIODOGRAMS; k++) {
        for (i = 1; i < LINEKEY_SIZE; i++) {
            for (m = (int)mel_bins[i-1]; m < (int)mel_bins[i]; m++) {
                after_mel[k][i] += periodograms[k][m];
            }
        }
    }
    
    /* average on all mel results -> At the very end + Decibel */
    memset(psd, 0, sizeof(psd));
    for (i = 0; i < LINEKEY_SIZE; i++) {
        psd[i] = 0;
        for (k = 0; k < NUM_OF_PERIODOGRAMS; k++) {
            psd[i] += after_mel[k][i];
        }
        //psd[i] = 10 * log10(psd[i]);
    }
    
    /* Dynamic threashold */
    memset(linekey_bits, 0, sizeof(linekey_bits));
    if (threshold_function(psd, &a, &b) < 0) {
        return -1;
    }
    
    for (i = 0; i < LINEKEY_SIZE; i++) {
        /* convert psd[i] to dB */
        psd[i] = 10 * log10(psd[i]);
        /* check threshold value */
        cur_threshold = (a + (b * (i+1)))*(10.0 / log(10)) + THRESHOLD_BIAS;
        linekey_bits[i] = psd[i] > cur_threshold;
        //printf("PSD at %d: value: %lf :: threshold: %lf -> bit=%d\n", i, psd[i], cur_threshold, linekey_bits[i]);
    }
    
    /* Final result */
    *out_linekey = 0;
    //printf("LineKey: ");
    for (i = 0; i < LINEKEY_SIZE; i++) {
        // printf("%d", linekey_bits[i]);
        *out_linekey = (*out_linekey << 1) | (linekey_bits[i] & 0x01);
    }
    

  out: 
    /* Cleanup */
    if (periodograms) {
        for (i = 0; i < NUM_OF_PERIODOGRAMS; i++) {
            if (periodograms[i]) free(periodograms[i]);
        }
        free(periodograms);
    }
    if (fft_plan) fftw_destroy_plan(fft_plan);
    if (out) fftw_free(out);

    return return_val;
}

int anlyze_new_data(const double* audio_data, const int data_size, const wav_header_t header, const int song_idx) {
    
    /* General variables */
    int return_val = 0;
    int i = 0;

    /* Window data - Calculations */
    int ms_size = (header.sample_rate * header.num_channels) / 1000; // represent the number of samples in a single ms
    int sigma_size = SIGMA_WINDOW_LENGTH * ms_size;
    int sub_size = SUB_WINDOW_LENGTH * ms_size;
    int step_size = STEP_LENGTH * ms_size;

    double cur_window[sigma_size];
    uint64_t cur_linekey = 0;
    
    if (!audio_data) {
        return_val = -1;
        goto out;
    }

    /* Window Calculations */
    for (i = 0; i <= data_size-step_size; i += step_size) {
        /* set current window */
        memcpy(cur_window, audio_data + i, sigma_size);

        /* Generate LineKey */
        new_linekey_generation(cur_window, sigma_size, sub_size, header.sample_rate, &cur_linekey);
        
        /* Save LineKey in DataBase */
        if (new_linekey_entry(cur_linekey, song_idx, i) < 0) {
            fprintf(stderr, "Error: Failed ot enter linekey for song %d at position %d.\n", song_idx, i);
            return_val = -1;
            goto out;
        }
    }
    
  out:
    return return_val;
}

int read_and_convert(char *path_to_file, double **sample_output, wav_header_t *out_header, long *size) {
    unsigned char* data_buffer= NULL;
    int read_num;

    /* Read raw file */
    read_num = read_file(path_to_file, &data_buffer, size, out_header);
    if (!data_buffer || read_num != 0) {
        fprintf(stderr, "Error getting data from file.\n");
        if (data_buffer) free(data_buffer);
        return -1;
    }

    /* Convert to samples */
    convert_to_samples(data_buffer, sample_output, size, out_header->bits_per_sample);
    /* Free data -> no need for the raw data anymore */
    if (data_buffer) free(data_buffer); // Free space on RAM.
    return 0;
}

int record_audio(char *output_path) {
    int status = 0;
    char command[256];
    
    memset(command, 0, sizeof(command));

    snprintf(command, sizeof(command), "arecord -d %d -f S16_LE -r 44100 -c 1 %s", RECORD_DURATION, output_path);

    printf("Recording...\n");

    status = system(command);

    if (!status) {
        printf("Recording saved to %s\n", output_path);
        return 0;
    } else {
        fprintf(stderr, "Error: Recording failed.\n");
        return -1;
    }
}

uint64_t key_dist(const uint64_t key_1, const uint64_t key_2) {
    uint64_t xor = 0;
    uint64_t dist = 0;
    int i = 0;

    /* XOR to get different bits (the distance) */
    xor = key_1 ^ key_2;

    /* Count how many 1s are there - the actual distance */
    for (i = 0; i < 64; i++)
        if ((xor >> i) & 1)
            dist++;

    return dist;
}

uint64_t get_min_dist(const uint64_t linekey) {
    FILE *f = NULL;
    uint64_t cur_dist = 0;
    uint64_t min_dist = -1; // Max value for uint64_t.
    linekey_t cur_key;

    f = fopen(PATH_TO_L, "rb");
    if (!f) return -1;

    while (fread(&cur_key, sizeof(cur_key), 1, f) == 1) {
        cur_dist = key_dist(linekey, cur_key.value);
        if (cur_dist < min_dist) {
            min_dist = cur_dist;
        }
    }

    if (f) fclose(f);
    f = NULL;
    return min_dist;
}

void add_new_element(dynamic_list_t *list, linekey_t new_element) {
    if (!list) return;
    void *temp = NULL;
    if (list->count == list->capacity) {
        /* Double the length */
        list->capacity = (list->capacity == 0) ? 4 : list->capacity * 2;
        temp = realloc(list->data, list->capacity * sizeof(linekey_t));
        if (!temp) return;
        list->data = temp;
    }
    list->data[list->count++] = new_element;
}

int get_closest_list(const uint64_t linekey, dynamic_list_t *key_list) {
    uint64_t min_dist = 0;
    FILE *f = NULL;
    linekey_t cur_key;

    f = fopen(PATH_TO_L, "rb");
    if (!f) return -1;

    min_dist = get_min_dist(linekey);
    while (fread(&cur_key, sizeof(cur_key), 1, f) == 1) {
        if (key_dist(linekey, cur_key.value) <= min_dist) {
            add_new_element(key_list, cur_key);
        }
    }

    if (f) fclose(f);
    f = NULL;
    return 0;
}

int get_majority(const int *arr, const int len) {
    int i = 0;
    int max = 0;
    int *backets = NULL;
    int max_count = 0;
    int majority = -1;

    /* Find max for range of backets */
    for (i = 0; i < len; i++) {
        if (max < arr[i])
            max = arr[i];
    }

    /* allocate space */
    backets = calloc(max+1, sizeof(int));

    /* Find majority */
    for (i = 0; i < len; i++) {
        backets[arr[i]]++;
        if (backets[arr[i]] > max_count) {
            max_count = backets[arr[i]];
            majority = arr[i];
        }
    }

    /* Cleanup */
    free(backets);

    return majority;
}

int recognize_recording(char *output_path) {
    /* General Variables */
    int return_val = 0;
    int i = 0;
    int j = 0;
    int k = 0;

    /* Recorded Audio Variables */
    double* audio_data = NULL;
    wav_header_t header;
    long size;

    /* Analyze Variables */
    int num_of_keys = (RECORD_DURATION*1000 - SIGMA_WINDOW_LENGTH) / STEP_LENGTH + 1;
    uint64_t *rec_keys = malloc(num_of_keys * sizeof(uint64_t));
    int ms_size = 0;
    int sigma_size = 0;
    int sub_size = 0;
    int step_size = 0;

    /* Recognize Variables */
    dynamic_list_t *dist_mat[num_of_keys]; // each element is a list of LineKeys
    linekey_info_t *cur_list = NULL;
    int cur_len = 0;
    int majority_list[num_of_keys];
    /* Majority of L_i_j */
    int *id_list = NULL;
    int maj_id;
    /* Majority of mat[i] */
    int *id_list_i = NULL;

    if (record_audio(output_path) < 0) {
        fprintf(stderr, "Error: Could not record audio.\n");
        free(rec_keys);
        return -1;
    }

    /* Read file */
    read_and_convert(output_path, &audio_data, &header, &size);

    /* Window data - Calculations */
    ms_size = (header.sample_rate * header.num_channels) / 1000; // represent the number of samples in a single ms
    sigma_size = SIGMA_WINDOW_LENGTH * ms_size;
    sub_size = SUB_WINDOW_LENGTH * ms_size;
    step_size = STEP_LENGTH * ms_size;

    double cur_window[sigma_size]; // I want it on the stack

    /* Get all keys */
    if (!audio_data) {
        return_val = -1;
        goto out;
    }

    /* Check recording size */
    if (size / ms_size >= RECORD_DURATION*1000 + 100) {
        fprintf(stderr, "Error: Recording is longer than 5 seconds at %ld.\n", size / ms_size);
        return_val = -1;
        goto out;
    }

    /* Init distance matrix */
    for (i = 0; i < num_of_keys; i++) {
        dist_mat[i] = malloc(sizeof(dynamic_list_t));
        dist_mat[i]->count = 0;
        dist_mat[i]->capacity = 0;
        dist_mat[i]->data = NULL;
    }

    /* Window Calculations */
    for (i = 0; i <= size-step_size && j < num_of_keys; i += step_size) {
        /* set current window */
        memcpy(cur_window, audio_data + i, sigma_size);

        /* Generate LineKey */
        new_linekey_generation(cur_window, sigma_size, sub_size, header.sample_rate, &rec_keys[j]);

        /* Get distant list */
        get_closest_list(rec_keys[j], dist_mat[j]);
        j++;
    }

    /* Search song in database */
    for (i = 0; i < num_of_keys; i++) {
        /* dist_mat[i] -> fynamic list of linekeys */
        id_list_i = calloc(dist_mat[i]->count, sizeof(int));
        for (j = 0; j < dist_mat[i]->count; j++) {
            /* Get list of info */
            read_linekey_list(dist_mat[i]->data[j].value, &cur_list, &cur_len);
            /* Get int list of IDs */
            if (cur_len > 0) id_list = calloc(cur_len, sizeof(int));
            for (k = 0; k < cur_len; k++) {
                id_list[k] = cur_list[k].id;
            }
            /* Call function that finds majority */
            maj_id = get_majority(id_list, cur_len);
            id_list_i[j] = maj_id;   
            /* Cleanup */
            if (cur_list) free(cur_list);
            if (id_list) free(id_list);
            cur_len = 0;
        }
        /* Get i majority */
        maj_id = get_majority(id_list_i, dist_mat[i]->count);
        majority_list[i] = maj_id;
        /* Cleanup */
        free(id_list_i);
    }

    /* Get song index */
    return_val = get_majority(majority_list, num_of_keys);
    
  out:
    for (i = 0; i < num_of_keys; i++) {
        if (dist_mat[i]) {
            free(dist_mat[i]->data);
            free(dist_mat[i]);
        }
    }
    free(rec_keys);
    if (audio_data) free(audio_data);
    return return_val;
}
