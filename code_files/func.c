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

// int analyze_audio(char* data, long size, wav_header_t header) {
//     if (!data || size==0) {
//         fprintf(stderr, "Error: data not received.\n");
//         return 1;
//     }
//     int ms_size = (header.sample_rate * header.num_channels * header.bits_per_sample) / 1000;
//     int sigma_size = SIGMA_WINDOW_LENGTH * ms_size;
//     int sub_size = SUB_WINDOW_LENGTH * ms_size;
//     int step_size = STEP_LENGTH * ms_size;

//     char* cur_window = NULL;
//     int temp_linekey;

//     int data_offset = 0;
//     while (data_offset < size) {
//         cur_window = data+data_offset;
//         // analyze current window -> should generate a linekey!
//         analyze_window(cur_window, sigma_size, sub_size, header.bits_per_sample, &temp_linekey);
//     }

//     return 0;
// }

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
    int i=0;

    if (!sub_window || !sub_size) {
        return 1;
    }

    for (i=0; i < sub_size; i++) {
        hamming_value = 0.54 - 0.46 * cos(2 * M_PI * i / (sub_size-1));
        sub_window[i] *= hamming_value;
    }

    return 0;
}

int create_periodogram(double* sub_window, int sub_size) {
    // receiving a window after hamming calculation
    // now it needs to got through FFT
    // output is a periodogram -> power of FFT result
    return 0;
}

int generate_linekey(double* window_data, int sigma_size, int sub_size, linekey_t* out_linekey) {
    //variable
    double** periodograms = NULL;
    double* sub_window = NULL;
    int step_size=0;
    int i=0;

    if (!window_data || !sub_size || !sigma_size) {
        return 1;
    }
    // sub_size is the number of "points" (samples) in the FFT calculation (or hamming for that matter)
    step_size = ((sigma_size - sub_size) / NUM_OF_PERIODOGRAMS);
    periodograms = (double**)calloc(NUM_OF_PERIODOGRAMS, sizeof(double*));
    sub_window = window_data; // start position of the first sub window

    // Creation of Periodograms
    for (i=0; i<NUM_OF_PERIODOGRAMS; i++) {
        sub_window += (i*step_size);
        hamming(sub_window, sub_size);
    }

    // cleanup
    for (i=0; i<NUM_OF_PERIODOGRAMS; i++) {
        free(periodograms[i]);
    }
    free(periodograms);
    return 0;
}

int analyze_data(double* audio_data, long data_size, wav_header_t header, int song_idx, char** unique_linekeys) {
    int ms_size = 0;
    int sigma_size = 0;
    int sub_size = 0;
    int step_size = 0;

    if (!audio_data || !data_size) {
        return 1;
    }

    ms_size = (header.sample_rate * header.num_channels) / 1000; // represent the number of samples in a single ms
    sigma_size = SIGMA_WINDOW_LENGTH * ms_size;
    sub_size = SUB_WINDOW_LENGTH * ms_size;
    step_size = STEP_LENGTH * ms_size;

    return 0;
}

// int hamming(char* sub_window, double* cur_hamming, int sub_size, int bits_per_sample) {
//     if (!sub_window || sub_size==0) {
//         return 1;
//     }
//     // Here we calculate the hamming function for each sample
//     int bytes_per_sample = bits_per_sample / 8;
//     for (int i=0; i < sub_size; i+= bytes_per_sample) {
//         double humming_value = 0.54 - 0.46 * cos(2 * M_PI * ((double)i/bytes_per_sample) / (NUM_POINTS - 1));
//         cur_hamming[i/bytes_per_sample] = *(sub_window+i) * humming_value;
//     }
//     //double hamming_value = 0.54 - 0.46 * cos(2 * M_PI * i / (NUM_POINTS - 1));
//     return 0;
// }

// int analyze_window(char* window_data, int sigma_size, int sub_size, int bits_per_sample, int* linekey) {
//     if (!window_data) {
//         fprintf(stderr, "Error: Window not received.\n");
//         return 1;
//     }
    
//     // define length and numbers
//     int periodogram_num = 16;
//     int step_size = ((sigma_size - sub_size) / periodogram_num);

//     // Loop on sub windows
//     // Each run of the loop computes the periodogram for sub window
//     for (int i=0; i < sigma_size; i+=step_size) {
//         char* cur_sub_window = window_data + i;
//         // The size of this array is the number of samples in each sub window
//         // This number is calculated to be the size of sub_window (in bytes)
//         // divided by the number of bytes in each sample (bits_per_sample / 8)
//         double cur_hamming[sub_size/(bits_per_sample/8)];
//         hamming(cur_sub_window, cur_hamming, sub_size, bits_per_sample);
//     }

//     return 0;
// }
