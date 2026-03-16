#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/main.h"

int main() {
    unsigned char* data_buffer= NULL;
    long size = 0;
    wav_header_t header;
    double* audio_data = NULL;
    /*
    linekey_t* unique_linekeys = NULL;
    int unique_len = 0;
    int song_idx = 0;
    int i = 0;
    uint64_t cur_linekey = 0;
    */

    int read_num = read_file(PATH_TO_FILE, &data_buffer, &size, &header);
    if (!data_buffer || read_num != 0) {
        fprintf(stderr, "Error getting data from file.\n");
        return 1;
    }

    convert_to_samples(data_buffer, &audio_data, &size, header.bits_per_sample);

    int ms_size = (header.sample_rate * header.num_channels) / 1000; // represent the number of samples in a single ms
    int sigma_size = SIGMA_WINDOW_LENGTH * ms_size;
    int sub_size = SUB_WINDOW_LENGTH * ms_size;
    printf("Size of ms: %d :: Size of window: %d || Size of sub window: %d\n", ms_size, sigma_size, sub_size);
    uint64_t linekey;
    new_linekey_generation(audio_data+(sigma_size*3), sigma_size, sub_size, header.sample_rate, &linekey);
    printf("LineKey of a single window: %lu\n", linekey);
    
    free(data_buffer);
    
    /*
    // prepare arrays
    unique_linekeys = (linekey_t*)calloc(500, sizeof(linekey_t));
    unique_len = 500;
    for (i = 0; i < unique_len; i++) {
        unique_linekeys[i].position_arr = (int**)calloc(20, sizeof(int*));
        unique_linekeys[i].pos_arr_len = 20;
        unique_linekeys[i].new_pos = 0;
    }

    analyze_data(audio_data, (const) size, header, song_idx, unique_linekeys, &unique_len);

    printf("Final length of unique array is: %d\n", unique_len);

    for (i = 0; i < unique_len; i++) {
        cur_linekey = convert_linekey_to_number(unique_linekeys[i].binary_values);
        printf("linekey[%d] = %lu\n", i, cur_linekey);
    }

    printf("Done generation of linekeys, please continue to making the database.\n");

    free(unique_linekeys);
    free(data_buffer);
    */
    free(audio_data);
    return 0;
}