#ifndef __FUNC_H__
#define __FUNC_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h> // For fixed-width integer types
#include <inttypes.h>
#include <math.h>
#include <fftw3.h>
#include "database.h"

/* Constants */
#define SIGMA_WINDOW_LENGTH 300
#define SUB_WINDOW_LENGTH   200
#define STEP_LENGTH         100
#define NUM_OF_PERIODOGRAMS 16
// Mel filter bank values
#define MIN_FREQ            0
#define MAX_FREQ            44100
// linekey values
#define LINEKEY_SIZE        64
#define THRESHOLD_BIAS      3
/* Recording */
#define RECORD_DURATION     5

/* Structures and Data Types */
typedef struct wav_header {
    // RIFF Chunk Descriptor
    char     riff_header[4];        // "RIFF"
    uint32_t chunk_size;            // File size in bytes - 8
    char     wave_header[4];        // "WAVE"

    // "fmt " Sub-chunk
    char     fmt_header[4];         // "fmt "
    uint32_t fmt_chunk_size;        // Size of the fmt chunk (usually 16)
    uint16_t audio_format;          // Audio format (1 for PCM)
    uint16_t num_channels;          // Number of channels (1 for mono, 2 for stereo)
    uint32_t sample_rate;           // Sampling frequency (e.g., 44100)
    uint32_t byte_rate;             // sample_rate * num_channels * bits_per_sample / 8
    uint16_t block_align;           // num_channels * bits_per_sample / 8
    uint16_t bits_per_sample;       // 8, 16, 24, etc.
} wav_header_t;

typedef struct dynamic_list {
    linekey_t *data;
    int count;
    int capacity; // How much we have alocated.
} dynamic_list_t;

/* Library Functions */
int read_file(const char* input_file, unsigned char** out, long* size, wav_header_t* header);
int read_and_convert(char *path_to_file, double **sample_output, wav_header_t *out_header, long *size);
int convert_to_samples(unsigned char* data_buffer, double** sample_data, long* size, int bits_per_sample);
int new_linekey_generation(double* window_data, int sigma_size, int sub_size, int sample_rate, uint64_t* out_linekey);
int anlyze_new_data(const double* audio_data, const int data_size, const wav_header_t header, const int song_idx);
int recognize_recording(char *output_path);

#endif