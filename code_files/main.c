#include "main.h"
#include "func.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main() {
    unsigned char* data_buffer= NULL;
    double* audio_data = NULL;
    linekey_t* unique_linekeys = NULL;
    long size = 0;
    wav_header_t header;

    int read_num = read_file(PATH_TO_FILE, &data_buffer, &size, &header);
    if (!data_buffer || read_num != 0) {
        fprintf(stderr, "Error getting data from file.\n");
        return 1;
    }

    convert_to_samples(data_buffer, &audio_data, &size, header.bits_per_sample);

    return 0;
}