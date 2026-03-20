#include "../include/database.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <threads.h>

int add_new_c_entry(linekey_info_t new_info, uint64_t start_offset) {
    /* NOTE: The function does not close the file! */
    int return_val = 0;
    uint32_t cur_pos = 0;
    uint32_t last_pos= 0;
    uint32_t end_pos = 0;
    linekey_info_t cur_info;
    FILE *f = NULL;
    memset(&cur_info, 0, sizeof(cur_info));

    f = fopen(PATH_TO_C, "rb+");
    if (!f) {
        /* File does not exist */
        f = fopen(PATH_TO_C, "wb+");
        if (!f) return -1;
    }

    /* Check if empty */
    fseek(f, 0, SEEK_END);
    if (!ftell(f) || !start_offset) {
        // if (!start_offset && ftell(f)) {
        //     /* Creating new Database -> clear old one. */
        //     fclose(f);
        //     f = NULL;
        //     f = fopen(PATH_TO_C, "wb+");
        //     if (!f) return -1;
        // }
        
        /* File is empty -> write new info at the start and exit */
        fseek(f, 0, SEEK_SET);
        fwrite(&new_info, sizeof(new_info), 1, f);
        fflush(f);
        return_val = 0;
        goto out;
    }
    
    /* Go to start offset in file */
    cur_pos = start_offset;
    fseek(f, cur_pos, SEEK_SET);

    while (fread(&cur_info, sizeof(cur_info), 1, f) == 1) {
        if (cur_info.next_offset == 0) {
            /* Found the last node of the list -> append new one */
            last_pos = ftell(f) - sizeof(linekey_info_t);
            fseek(f, 0, SEEK_END);
            end_pos = ftell(f);
            /* Write new node */
            fwrite(&new_info, sizeof(new_info), 1, f);
            fflush(f);
            /* Link the last node to the list */
            cur_info.next_offset = end_pos;
            fseek(f, last_pos, SEEK_SET);
            fwrite(&cur_info, sizeof(cur_info), 1, f);
            fflush(f);

            return_val = 0;
            goto out;
        }
        cur_pos = cur_info.next_offset;
        fseek(f, cur_pos, SEEK_SET);
    }
    /* Failed to add new info */
    return_val = -1;

  out:
    if (f) fclose(f);
    f = NULL;
    return return_val;
}

int new_linekey_entry(const uint64_t new_linekey, const int song_idx, const int position) {
    int return_val = 0;
    FILE *fl = NULL;
    FILE *fc = NULL;
    linekey_t cur_linekey;
    linekey_info_t info;
    memset(&cur_linekey, 0, sizeof(cur_linekey));
    memset(&info, 0, sizeof(info));

    /* Init new info */
    info.id = song_idx;
    info.p = position;
    info.next_offset = 0;

    /* Open file - if not exists create */
    fl = fopen(PATH_TO_L, "rb+");
    if (!fl) {
        /* File does not exist - open a new one */
        fl = fopen(PATH_TO_L, "wb+");
        /* new open failed! */
        if (!fl) {
            fprintf(stderr, "Error: Failed to open L file.\n");
            return_val = -1;
            goto out;
        }
        cur_linekey.start_offset = 0;
        cur_linekey.value = new_linekey;
        cur_linekey.count = 1;
        fwrite(&cur_linekey, sizeof(cur_linekey), 1, fl);
        fflush(fl);
        /* Write to C file */
        return_val = add_new_c_entry(info, cur_linekey.start_offset);
        if (return_val < 0)
            fprintf(stderr, "Error: Failed to insert data to C file.\n");
        goto out;
    } else {
        /* Find the linekey */
        while (fread(&cur_linekey, sizeof(cur_linekey), 1, fl) == 1) {
            if (cur_linekey.value == new_linekey) {
                cur_linekey.count++;
                /* Go back to start of cur_linekey and rewrite it */
                fseek(fl, -sizeof(cur_linekey), SEEK_CUR);
                fwrite(&cur_linekey, sizeof(cur_linekey), 1, fl);
                fflush(fl);
                /* Update the list in C */
                return_val = add_new_c_entry(info, cur_linekey.start_offset);
                if (return_val < 0)
                    fprintf(stderr, "Error: Failed to insert data to C file.\n");
                goto out;
            }
        }

        /* Not found the linekey in the database -> add new entry */
        fseek(fl, 0, SEEK_END);
        cur_linekey.value = new_linekey;
        cur_linekey.count = 1;
        /* Start offset is the end of C file */
        fc = fopen(PATH_TO_C, "ab+");
        if (!fc) {
            /* Does not exist -> offset if 0 */
            cur_linekey.start_offset = 0;
        } else {
            fseek(fc, 0, SEEK_END);
            cur_linekey.start_offset = ftell(fc);
            /* Write new entry to C file */
            fwrite(&info, sizeof(info), 1, fc);
            fflush(fc);
            fclose(fc);
            fc = NULL;
        }
        fwrite(&cur_linekey, sizeof(cur_linekey), 1, fl);
        fflush(fl);
        goto out;
    }

  out:
    if (fl) fclose(fl);
    fl = NULL;
    if (fc) fclose(fc);
    fc = NULL;
    return return_val;
}

void print_linekey(linekey_t linekey) {
    printf("Value: 0x%016" PRIx64 " :: Offset: %lu :: Count: %lu\n", linekey.value, linekey.start_offset, linekey.count);
}

void print_info(linekey_info_t info, int cur_loc) {
    printf("ID: %d :: Position: %d :: Next: %d :: Location in file: %d\n", info.id, info.p, info.next_offset, cur_loc);
}

void print_metadata(meta_data_t data) {
    printf("name: %s\nalbum: %s\nartist: %s\n\n", data.name, data.album, data.artist);
}

int print_data(const char *name) {
    int return_val = 0;
    FILE *f = NULL;
    int type = 0; /* 0 -> linekey || 1 -> info */
    linekey_t cur_linekey = {0};
    linekey_info_t cur_info = {0};
    meta_data_t cur_data = {0};

    if (!name) return -1;
    
    /* Open file */
    if (!strncmp(name, L_NAME, sizeof(L_NAME))) {
        /* File is L */
        f = fopen(PATH_TO_L, "rb");
        if (!f) {
            return_val = -1;
            goto out;
        }
    } else if (!strncmp(name, C_NAME, sizeof(C_NAME))) {
        f = fopen(PATH_TO_C, "rb");
        if (!f) {
            return_val = -1;
            goto out;
        }
        type = 1;
    } else if (!strncmp(name, M_NAME, sizeof(M_NAME))) {
        f = fopen(PATH_TO_M, "rb");
        if (!f) {
            return_val = -1;
            goto out;
        }
        type = 2;
    } else {
        return_val = -1;
        goto out;
    }

    /* Read File -> currently print */
    if (type == 1) {
        while (fread(&cur_info, sizeof(cur_info), 1, f) == 1)
            print_info(cur_info, ftell(f) - sizeof(cur_info));
    } else if (type == 0){
        while (fread(&cur_linekey, sizeof(cur_linekey), 1, f) == 1)
            print_linekey(cur_linekey);
    } else {
        while (fread(&cur_data, sizeof(cur_data), 1, f) == 1)
            print_metadata(cur_data);
    }

  out:
    if (f) fclose(f);
    f = NULL;
    return return_val;
}

int reset_database(void) {
    FILE *f = NULL;

    f = fopen(PATH_TO_L, "wb");
    if (!f) return -1;
    fclose(f);
    f = NULL;
    f = fopen(PATH_TO_C, "wb");
    if (!f) return -1;
    fclose(f);
    f = NULL;
    f = fopen(PATH_TO_M, "wb");
    if (!f) return -1;
    fclose(f);
    f = NULL;
    return 0;
}

int get_linekey(const uint64_t linekey, linekey_t *out) {
    FILE *f = NULL;
    linekey_t cur_key;

    f= fopen(PATH_TO_L, "rb");
    if (!f) return -1;

    while (fread(&cur_key, sizeof(cur_key), 1, f) == 1) {
        if (cur_key.value == linekey) {
            *out = cur_key;
            return 0;
        }
    }
    return -1;
}

int read_linekey_list(const uint64_t linekey, linekey_info_t **out_list, int *out_len) {
    /* Get the list:
     * 1. Search LineKey in L file.
     * 2. Allocate space according to count variable.
     * 3. Go to start offset in C and start getting the list.
     */
    
    int return_val = 0;
    uint64_t i = 0;
    FILE *f = NULL;
    linekey_t key_obj;
    uint32_t next_addr = 0;

    return_val = get_linekey(linekey, &key_obj);
    if (return_val < 0) goto out;

    if (!out_len) return -1;
    *out_len = key_obj.count;
    if (!key_obj.count) return 0;
    *out_list = (linekey_info_t *) calloc(key_obj.count, sizeof(linekey_info_t));

    f = fopen(PATH_TO_C, "rb");
    if (!f) return -1;
    next_addr = key_obj.start_offset;
    for (i = 0; i < key_obj.count; i++) {
        fseek(f, next_addr, SEEK_SET);
        if (fread(&((*out_list)[i]), sizeof(linekey_info_t), 1, f) != 1) {
            return_val = -1;
            goto out;
        }
        next_addr = (*out_list)[i].next_offset;
        if (!next_addr && i < key_obj.count - 1) {
            return_val = -1;
            goto out;
        }
    }

  out:
    if (f) fclose(f);
    f = NULL;
    return return_val;
}

int insert_new_metadata(const meta_data_t *new_data) {
    FILE *f = NULL;

    f = fopen(PATH_TO_M, "ab+");
    if (!f) return -1;

    fwrite(new_data, sizeof(meta_data_t), 1, f);

    if (f) fclose(f);
    f = NULL;

    return 0;
}

int retrieve_metadata(int idx, meta_data_t * data) {
    FILE *f = NULL;

    f = fopen(PATH_TO_M, "rb");
    if (!f) return -1;

    fseek(f, sizeof(meta_data_t) * (idx-1), SEEK_SET);

    if (fread(data, sizeof(*data), 1, f) != 1)
        return -1;

    if (f) fclose(f);
    f = NULL;

    return 0;
}
