#include "../include/main.h"
#include <string.h>

int main() {
    
    // long size = 0;
    // wav_header_t header;
    // double* audio_data = NULL;
    /* Reset Database -> Remove later! */
    // reset_database();

    // /* Read raw file 2*/
    // read_and_convert(PATH_TO_FILE_2, &audio_data, &header, &size);

    // anlyze_new_data(audio_data, size, header, 2);

    /* Read database */
    // printf("L file:\n");
    // print_data(L_NAME);
    // printf("C file:\n");
    // print_data(C_NAME);

    // meta_data_t data = {0};
    // strcpy(data.name, "Love is a laser quest");
    // strcpy(data.album, "Suck it and see");
    // strcpy(data.artist, "Arctic Monkeys");
    // insert_new_metadata(&data);
    // memset(&data, 0, sizeof(data));

    // strcpy(data.name, "גשם חזק");
    // strcpy(data.album, "יחסים פתוחים");
    // strcpy(data.artist, "Monica Sex");
    // insert_new_metadata(&data);
    // memset(&data, 0, sizeof(data));
    
    // strcpy(data.name, "Piledriver walz");
    // strcpy(data.album, "Suck it and see");
    // strcpy(data.artist, "Arctic Monkeys");
    // insert_new_metadata(&data);    
    // memset(&data, 0, sizeof(data));
    
    // strcpy(data.album, "Whatever peaple say I am that's what I'm not");
    // strcpy(data.artist, "Arctic Monkeys");

    // /* Expand Database loop */
    // int lengths[NEW_SONGS_NUM] = {
    //     222,
    //     173,
    //     177,
    //     141,
    //     130,
    //     173,
    //     134,
    //     143,
    //     175,
    //     268,
    //     202,
    //     193,
    //     331
    // };
    // char names[NEW_SONGS_NUM][128] = {
    //     "The View From The Afternoon",
    //     "I Bet You Look Good On The Dance Floor",
    //     "Fake Tales Of San Francisco",
    //     "Dancing Shoes",
    //     "You Probably Couldn't See For The Light But You Were Staring Strait At Me",
    //     "Still Take You Home",
    //     "Riot Van",
    //     "Red Light Indicates Doors Are Secured",
    //     "Mardy Bum",
    //     "Perhaps Vampires Is A Bit Strong But...",
    //     "When The Sun Goes Down",
    //     "From The Ritz To The Rubble",
    //     "A Certain Romance"
    // };

    // for (int i = 0; i < NEW_SONGS_NUM; i++) {
    //     memset(&(data.name), 0, sizeof(data.name));
    //     strcpy(data.name, names[i]);
    //     insert_new_metadata(&data);
    //     /* Rewcord audio */
    //     record_audio(PATH_TO_FILE_1, lengths[i]);
    //     /* Read raw file 1*/
    //     read_and_convert(PATH_TO_FILE_1, &audio_data, &header, &size);

    //     anlyze_new_data(audio_data, size, header, i+4);

    //     if (audio_data) free(audio_data);
    //     audio_data = NULL;
    // }

    // printf("Meta Data File:\n");
    // print_data(M_NAME);

    /* Recognize a song */
    int id = recognize_recording(PATH_TO_TEMP_RECORDING);
    printf("Found song with idx %d.\n", id);
    meta_data_t data = {0};
    retrieve_metadata(id, &data);
    printf("Name: %s\nAlbum: %s\nArtist: %s\n", data.name, data.album, data.artist);

    /* Cleanup */
    // if (audio_data) free(audio_data);
    return 0;
}