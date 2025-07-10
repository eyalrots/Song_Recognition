#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include "func.h"

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
    int sample = 0;
    for (i=bytes_pet_sample-1; i>=0; i--) {
        sample = (sample << 8) | samples[i];
    }
    return (double)sample;
}

int convert_to_samples(unsigned char* data_buffer, double** sample_data, long* size, int bits_per_sample) {
    int i=0;
    int bytes_per_sample = bits_per_sample / 8;

    if (!data_buffer) {
        return 1;
    }

    *size = *size / bytes_per_sample;
    *sample_data = (double*)calloc(*size, sizeof(double));
    for (i=0; i < *size; i++) {
        (*sample_data)[i] = get_sample(&data_buffer[i*bytes_per_sample], bytes_per_sample);
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

void copy_arr(double* original_arr, double* output_arr, int size) {
    int i = 0;
    for (; i<size; i++) {
        output_arr[i] = original_arr[i];
    }
}

int create_periodogram(double* sub_window, int sub_size, int per_size, double* in, fftw_complex* out, fftw_plan p, double* periodogram) {
    // receiving a window after hamming calculation
    // now it needs to got through FFT
    // output is a periodogram -> power of FFT result
    int i = 0;
    if (!sub_window || !in || !out || !periodogram) {
        return 1;
    }
    // prepare input data for fft
    copy_arr(sub_window, in, sub_size);
    fftw_execute(p);
    for (i = 0; i < per_size; i++) {
        periodogram[i] = (out[i][0] * out[i][0] + out[i][1] * out[i][1]);
    }
    return 0;
}

int fft_to_mel(double* fft_values, double* mel_points, double delta_f) {
    double min_mel = 0;
    double max_mel = 0;
    int jump = 0;
    int i = 0;

    if (!fft_values || !mel_points) {
        return 1;
    }
    // map Hz to Mel
    min_mel = 2595 * log10(1 + ((double)MIN_FREQ / 700));
    max_mel = 2595 * log10(1 + ((double)MAX_FREQ / 700));

    // create points on mel scale spaced evenly
    jump = (max_mel - min_mel) / LINEKEY_SIZE;
    mel_points[0] = min_mel;
    for (i = 1; i < LINEKEY_SIZE+1; i++) {
        mel_points[i] = mel_points[i-1] + jump;
    }

    // check if space of points on mel scale was correct
    if (mel_points[LINEKEY_SIZE] <= max_mel-jump) {
        fprintf(stderr, "Error: wrong input to mel scale - end_of_scale != max_mel.\nActual value is: %lf while max is: %lf and the jump is: %d", mel_points[LINEKEY_SIZE], max_mel, jump);
        return 1;
    }

    // convert back to Hz
    for (i = 0; i < LINEKEY_SIZE+1; i++) {
        mel_points[i] = 700 * pow(10,((double)mel_points[i] / 2595) - 1);
    }

    // round to FFT values
    for (i = 0; i < LINEKEY_SIZE+1; i++) {
        mel_points[i] = (int)floor(mel_points[i] / delta_f);
    }

    return 0;
}

int pass_periodogram_through_mel_filter(double* periodogram, double* result, double start, double end) {
    int m = 0;

    if (!periodogram  || !result) {
        return 1;
    }

    *result = 0;
    for (m = start; m < end; m++) {
        *result += periodogram[m];
    }

    return 0;
}

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

int threshold_function(double* psd, double* out_func) {
    double sqr_sum = 0;
    double sum_sqr = 0;
    double ln_sum_for_b_numerator = 0;
    double denominator = 0;
    int i = 0;
    long double a_numerator = 0;
    long double b_numerator = 0;
    double a = 0;
    double b = 0;
    
    if (!psd || !out_func) {
        return 1;
    }

    // calculate squared sum
    for (i = 1; i < LINEKEY_SIZE+1; i++) {
        sqr_sum += i-1;
    }
    sqr_sum *= sqr_sum; // squaring the sum

    // calculate sum of squared values
    for (int i = 1; i < LINEKEY_SIZE + 1; i++) {
        sum_sqr += (i-1) * (i-1);
    }

    // calculate numerator
    denominator = LINEKEY_SIZE * sum_sqr - sqr_sum;

    // calculate a numerator
    a_numerator = 0;
    for (i = 0; i <LINEKEY_SIZE; i++) {
        a_numerator += log(psd[i]) * (sum_sqr - i);
        printf("a(t_n) numerator += [ln(%lf)=%lf] * (%lf - %d) = %Lf\n", psd[i], log(psd[i]), sum_sqr, i, a_numerator);
    }
    printf("a(t_n) numerator = %Lf\n", a_numerator);

    // calculate ln sum for b numerator
    ln_sum_for_b_numerator = 0;
    for (i = 0; i < LINEKEY_SIZE; i++) {
        ln_sum_for_b_numerator += log(psd[i]);
    }

    // calculate b numerator
    b_numerator = 0;
    for (i = 0; i < LINEKEY_SIZE; i++) {
        b_numerator += i * (LINEKEY_SIZE * log(psd[i]) - ln_sum_for_b_numerator);
    }

    // calculate a(t_n) and b(t_n)
    a = a_numerator / denominator;
    b = b_numerator / denominator;
    printf("a(t_n) = %Lf / %lf = %lf\n", a_numerator, denominator, a);
    printf("b(t_n) = %Lf / %lf = %lf\n", b_numerator, denominator, b);

    // calculate out function
    for (i = 0; i < LINEKEY_SIZE; i++) {
        out_func[i] = (a * exp(b * i)) + THRESHOLD_BIAS; // eq. (6) in paper
    }

    return 0;
}

int generate_linekey(double* window_data, int sigma_size, int sub_size, int sample_rate, linekey_t* out_linekey) {
    //variable
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

    //prepare fftw plan
    per_size = (sub_size / 2 + 1);
    in = (double *) fftw_malloc(sizeof(double) * sub_size);
    out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * per_size);
    if (!in || !out) {
        fprintf(stderr, "Error: Memory allocation failed.\n");
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
        hamming(sub_window, sub_size);
        create_periodogram(sub_window, sub_size, per_size, in, out, p, periodograms[i]); // eq. (2) in paper
    }

    // calculate FFT values that are present in the PSD
    delta_f = (double)sample_rate / sub_size;
    fft_values = (double*) calloc(per_size, sizeof(double));
    for (i = 0; i < per_size; i++) {
        fft_values[i] = i * delta_f;
    }

    // Mel filter bank calculation
    mel_points = (double *) calloc(LINEKEY_SIZE+1, sizeof(double));
    mel_filter_results = (double**) calloc(NUM_OF_PERIODOGRAMS, sizeof(double*));
    fft_to_mel(fft_values, mel_points, delta_f);
    // calculate mel results for each periodogram
    for (i = 0; i < NUM_OF_PERIODOGRAMS; i++) {
        mel_filter_results[i] = (double*) calloc(LINEKEY_SIZE, sizeof(double));
        for (j = 1; j < LINEKEY_SIZE+1; j++) {
            pass_periodogram_through_mel_filter(periodograms[i], &mel_filter_results[i][j-1], mel_points[j-1], mel_points[j]); // eq. (3) in paper
        }
    }

    // calculate PSD by summing up the mel filters
    psd = (double*) calloc(LINEKEY_SIZE, sizeof(double));
    calculate_psd(mel_filter_results, psd);

    // Adaptive threshold function
    threshold = (double*) calloc(LINEKEY_SIZE, sizeof(double));
    threshold_function(psd, threshold);

    // generating the linekey
    printf("Generating linekey.\n");
    for (i = 0; i < LINEKEY_SIZE; i++) {
        out_linekey->binary_values = (int*)calloc(LINEKEY_SIZE, sizeof(int));
        printf("linekey[%d] = (%lf > %lf)\n", i, psd[i], threshold[i]);
        out_linekey->binary_values[i] = (psd[i] > threshold[i]);
    }
    printf("Done generating linekey.\n");

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

int analyze_data(double* audio_data, long data_size, wav_header_t header, int song_idx, linekey_t* unique_linekeys, int* unique_len) {
    int ms_size = 0;
    int sigma_size = 0;
    int sub_size = 0;
    int step_size = 0;
    double* cur_window = NULL;
    int return_val = 0;
    int i = 0;
    int cur_linekey = 0;

    if (!audio_data || !data_size || !unique_linekeys) {
        fprintf(stderr, "Error: one of the input values is empty.\n");
        return_val = 1;
        goto out;
    }

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
        // expand unique linekey array if full
        if (cur_linekey >= *unique_len) {
            unique_linekeys = (linekey_t*)realloc(unique_linekeys, (*unique_len+50) * sizeof(linekey_t));
            *unique_len += 50;
        }
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
