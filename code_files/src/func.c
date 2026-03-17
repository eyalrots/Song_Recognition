#include "../include/func.h"

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

/*---Deprecated---*/
double triangular_filter(double f, double start, double peak, double end) {
    if (f <= start || f >= end) {
        return 0;
    }

    if (f >= start && f < peak) {
        return (f - start) / (peak - start);
    } else if (f >= peak && f <end) {
        return (end - f) / (end - peak);
    }

    return 0;
}

/*---Deprecated---*/
int pass_periodogram_through_triangular_filter(double* periodogram, double* result, int size, double start, double peak, double end) {
    int i = 0;
    
    if (!periodogram || !result) {
        return 1;
    }
    *result = 0;
    for (i = 0; i < size; i++) {
        *result += triangular_filter(periodogram[i], start, peak, end);
    }

    return 0;
}

/*---Deprecated---*/
int pass_periodogram_through_mel_filter(double* periodogram, double* result, double start, double end) {
    int m = 0;

    if (!periodogram  || !result) {
        return 1;
    }

    *result = 0;
    for (m = start; m < end; m++) {
        *result += periodogram[m];
    }

    //*result = *result / (end - start);

    return 0;
}

/*---Deprecated---*/
int calculate_psd(double** mel_filter_results, double* out_psd) {
    int i = 0;
    int j = 0;
    double per_sum = 0;

    if (!mel_filter_results || !out_psd) {
        return 1;
    }

    for (i = 0; i < LINEKEY_SIZE; i++) {
        per_sum = 0;
        for (j = 0; j < NUM_OF_PERIODOGRAMS; j++) {
            per_sum += mel_filter_results[j][i];
        }
        out_psd[i] = 10 * log10(per_sum);
        out_psd[i] = isinf(out_psd[i]) ? 1 : out_psd[i];
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
    printf("a: %lf :: b: %lf\n", a, b);
    for (i = 0; i < LINEKEY_SIZE; i++) {
        /* convert psd[i] to dB */
        psd[i] = 10 * log10(psd[i]);
        /* check threshold value */
        cur_threshold = (a + (b * (i+1)))*(10.0 / log(10)) + THRESHOLD_BIAS;
        linekey_bits[i] = psd[i] > cur_threshold;
        printf("PSD at %d: value: %lf :: threshold: %lf -> bit=%d\n", i, psd[i], cur_threshold, linekey_bits[i]);
    }
    printf("LineKey is created.\n");
    /* Final result */
    *out_linekey = 0;
    printf("LineKey: ");
    for (i = 0; i < LINEKEY_SIZE; i++) {
        printf("%d", linekey_bits[i]);
        *out_linekey = (*out_linekey << 1) | (linekey_bits[i] & 0x01);
    }
    printf("\n");

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

/*---Deprecated---*/
int generate_linekey(double* window_data, int sigma_size, int sub_size, int sample_rate, 
                        linekey_t* out_linekey) {
    // Variables
    double** periodograms = NULL;
    double* psd = NULL;
    double* sub_window = NULL;
    double* fft_values = NULL;
    double* mel_points = NULL;
    double* threshold = NULL;
    double** mel_filter_results = NULL; // total of LINEKEY_SIZE array each the result of a single mel filter
    double delta_f = 0;
    int step_size = 0;
    int per_size = 0;
    int i = 0;
    int j = 0;
    int return_val = 0;
    // FFT plan variables
    double *in = NULL;
    fftw_complex *out = NULL;
    fftw_plan p = NULL;

    if (!window_data || !sub_size || !sigma_size || !out_linekey) {
        return_val = 1;
        goto out;
    }

    // Prepare fftw plan
    per_size = (sub_size / 2 + 1);
    in = (double *) fftw_malloc(sizeof(double) * sub_size);
    out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * per_size);
    if (!in || !out) {
        fprintf(stderr, "Error: Memory allocation in FFTW failed.\n");
        return_val = 1;
        goto out;
    }
    // fftw plan creation
    p = fftw_plan_dft_r2c_1d(sub_size, in, out, FFTW_MEASURE);

    // sub_size is the number of "points" (samples) in the FFT calculation (or hamming for that matter)
    step_size = ((sigma_size - sub_size) / NUM_OF_PERIODOGRAMS);
    periodograms = (double**)calloc(NUM_OF_PERIODOGRAMS, sizeof(double*));
    sub_window = window_data; // start position of the first sub window

    // Creation of Periodograms
    for (i=0; i<NUM_OF_PERIODOGRAMS; i++) {
        sub_window += (i*step_size);
        periodograms[i] = (double *) calloc(per_size, sizeof(double));
        /* This is not good! it chabges the original samples!!!! */
        hamming(sub_window, sub_size);
        //create_periodogram(sub_window, sub_size, per_size, in, out, p, periodograms[i]); // eq. (2) in paper
    }

    // calculate FFT values that are present in the PSD
    delta_f = (double)sample_rate / per_size;
    fft_values = (double*) calloc(per_size, sizeof(double));
    for (i = 0; i < per_size; i++) {
        fft_values[i] = i * delta_f;
    }
    printf("last frequency: %lf\n", fft_values[per_size-1]);

    // Mel filter bank calculation
    mel_points = (double *) calloc(LINEKEY_SIZE+1, sizeof(double));
    mel_filter_results = (double**) calloc(NUM_OF_PERIODOGRAMS, sizeof(double*));
    fft_to_mel(3, mel_points, delta_f); // MOVE TO OUTER FUNCTION!!!
    // calculate mel results for each periodogram
    for (i = 0; i < NUM_OF_PERIODOGRAMS; i++) {
        mel_filter_results[i] = (double*) calloc(LINEKEY_SIZE, sizeof(double));
        for (j = 1; j < LINEKEY_SIZE+1; j++) {
            pass_periodogram_through_mel_filter(periodograms[i], &mel_filter_results[i][j-1], mel_points[j-1], mel_points[j]); // eq. (3) in paper
            //pass_periodogram_through_triangular_filter(periodograms[i], &mel_filter_results[i][j-1], per_size, mel_points[j-1], mel_points[j], mel_points[j+1]);
        }
    }

    // calculate PSD by summing up the mel filters
    psd = (double*) calloc(LINEKEY_SIZE, sizeof(double));
    calculate_psd(mel_filter_results, psd);

    // Adaptive threshold function
    threshold = (double*) calloc(LINEKEY_SIZE, sizeof(double));
    //threshold_function(psd, threshold);

    // generating the linekey
    /*
    printf("Generating linekey.\n");
    for (i = 0; i < LINEKEY_SIZE; i++) {
        out_linekey->binary_values = (int*)calloc(LINEKEY_SIZE, sizeof(int));
        printf("linekey[%d] = (%lf > %lf)\n", i, psd[i], threshold[i]);
        out_linekey->binary_values[i] = (psd[i] > threshold[i]);
    }
    printf("Done generating linekey.\n");
    */

out:// cleanup
    for (i=0; i<NUM_OF_PERIODOGRAMS; i++) {
        free(periodograms[i]);
    }
    free(periodograms);
    for (i=0; i<NUM_OF_PERIODOGRAMS; i++) {
        free(mel_filter_results[i]);
    }
    free(mel_filter_results);
    free(psd);
    free(threshold);
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
    free(mel_points);
    return return_val;
}

int new_anlyze_data(const double* audio_data, const int data_size, const wav_header_t header, int song_idx) {
    
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
    for (i = 0; i < data_size; i += step_size) {
        /* set current window */
        memcpy(cur_window, audio_data+i, sigma_size);

        /* Generate LineKey */
        new_linekey_generation(cur_window, sigma_size, sub_size, header.sample_rate, &cur_linekey);
        printf("Generated LineKey: %016" PRIx64 "\n", cur_linekey);
        
        /* Save LineKey in DataBase */
    }
    
  out:
    return return_val;
}

int analyze_data(double* audio_data, const uint32_t data_size, wav_header_t header, 
                    int song_idx, linekey_t* unique_linekeys, int* unique_len) {
    int ms_size = 0;
    int sigma_size = 0;
    int sub_size = 0;
    int step_size = 0;
    double* cur_window = NULL;
    int return_val = 0;
    int i = 0;
    int cur_linekey = 0;

    if (!audio_data || !data_size) {
        fprintf(stderr, "Error: one of the input values is empty.\n");
        return_val = 1;
        goto out;
    }
    /* I don't think likey array should be here! */
    if (!unique_linekeys && *unique_len) {
        fprintf(stderr, "Error: unique linekey length != 0 but pointer is null.\n");
        return_val = 1;
        goto out;
    } else if (!unique_linekeys) {
        unique_linekeys = calloc(50, sizeof(linekey_t));
        *unique_len = 50;
    }

    ms_size = (header.sample_rate * header.num_channels) / 1000; // represent the number of samples in a single ms
    sigma_size = SIGMA_WINDOW_LENGTH * ms_size;
    sub_size = SUB_WINDOW_LENGTH * ms_size;
    step_size = STEP_LENGTH * ms_size;

    // go through data with sigma windows
    for (i = 0; i < data_size; i+=step_size) {
        cur_window = audio_data + i;
        // expand unique linekey array if full -> maybe not here?
        if (cur_linekey >= *unique_len) {
            unique_linekeys = (linekey_t*)realloc(unique_linekeys, (*unique_len+50) * sizeof(linekey_t));
            *unique_len += 50;
        }
        /*
        generate_linekey(cur_window, sigma_size, sub_size, header.sample_rate, &unique_linekeys[cur_linekey]);
        printf("Generated linekey number %d\n", cur_linekey);
        // insert current data to position array (song index and position in said song)
        // check if there is an error with the array
        if (!unique_linekeys[cur_linekey].position_arr && unique_linekeys[cur_linekey].pos_arr_len) {
            fprintf(stderr, "Error: position array length != 0 but pointer is null.\n");
            return_val = 1;
            goto out;
        } 
        // if array is empty (first song input) then create a new array
        else if (!unique_linekeys[cur_linekey].position_arr) {
            unique_linekeys[cur_linekey].position_arr = (int**)calloc(20, sizeof(int*));
            unique_linekeys[cur_linekey].pos_arr_len = 20;
        } 
        // expand array if full
        else if (unique_linekeys[cur_linekey].new_pos >= unique_linekeys[cur_linekey].pos_arr_len) {
            unique_linekeys[cur_linekey].position_arr = (int**)realloc(unique_linekeys[cur_linekey].position_arr, (unique_linekeys[cur_linekey].pos_arr_len + 20) * sizeof(int*));
            unique_linekeys[cur_linekey].pos_arr_len += 20;
        }
        // insert to array the data
        unique_linekeys[cur_linekey].position_arr[unique_linekeys[cur_linekey].new_pos] = (int*) calloc(2, sizeof(int));
        unique_linekeys[cur_linekey].position_arr[unique_linekeys[cur_linekey].new_pos][0] = song_idx;
        unique_linekeys[cur_linekey].position_arr[unique_linekeys[cur_linekey].new_pos][1] = song_idx;
        // advance indices of arrays
        cur_linekey++;
        unique_linekeys[cur_linekey].new_pos++;
        */
    }    

out:
    return return_val;
}

uint64_t convert_linekey_to_number(int* linekey) {
    uint64_t result = 0;
    int i = 0;

    for (i = 0; i < LINEKEY_SIZE; i++) {
        printf("bit[%d] = %d\n", i, linekey[i]);
        result += linekey[i] * pow(2, i);
    }

    return result;
}
