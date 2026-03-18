#include "../include/main.h"

int main() {
    long size = 0;
    wav_header_t header;
    double* audio_data = NULL;
    /* Reset Database -> Remove later! */
    reset_database();

    /* Read raw file 1*/
    read_and_convert(PATH_TO_FILE_1, &audio_data, &header, &size);

    anlyze_new_data(audio_data, size, header, 1);

    if (audio_data) free(audio_data);

    /* Read raw file 1 again */
    read_and_convert(PATH_TO_FILE_2, &audio_data, &header, &size);
    
    /* Insert new data to database */
    anlyze_new_data(audio_data, size, header, 2);

    /* Read database */
    printf("L file:\n");
    print_data(L_NAME);
    printf("C file:\n");
    print_data(C_NAME);

    /* Cleanup */
    if (audio_data) free(audio_data);
    return 0;
}