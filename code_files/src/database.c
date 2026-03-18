#include "../include/database.h"
#include <stdint.h>
#include <stdio.h>
#include <string.h>

int search_linekey(uint64_t linekey, int *lk_idx) {

    int return_val = 0;
    FILE *f = NULL;
    uint64_t cur_linekey = 0;

    /* Open file - if not exists create */
    f = fopen(PATH_TO_L, "rb+");
    if (!f) {
        return_val = -1;
        goto out;
    }



  out:
    if (f) fclose(f);
    f = NULL;
    return return_val;
}

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
        if (!start_offset && ftell(f)) {
            /* Creating new Database -> clear old one. */
            fclose(f);
            f = NULL;
            f = fopen(PATH_TO_C, "wb+");
            if (!f) return -1;
        }
        /* File is empty -> write new info at the start and exit */
        fseek(f, 0, SEEK_SET);
        fwrite(&new_info, sizeof(new_info), 1, f);
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
            /* Link the last node to the list */
            cur_info.next_offset = end_pos;
            fseek(f, last_pos, SEEK_SET);
            fwrite(&cur_info, sizeof(cur_info), 1, f);

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
            return_val = -1;
            goto out;
        }
        cur_linekey.start_offset = 0;
        cur_linekey.value = new_linekey;
        cur_linekey.count = 1;
        fwrite(&cur_linekey, sizeof(cur_linekey), 1, fl);
        /* Write to C file */
        return_val = add_new_c_entry(info, cur_linekey.start_offset);
        goto out;
    } else {
        /* Find the linekey */
        while (fread(&cur_linekey, sizeof(cur_linekey), 1, fl) == 1) {
            if (cur_linekey.value == new_linekey) {
                cur_linekey.count++;
                /* Go back to start of cur_linekey and rewrite it */
                fseek(fl, -sizeof(cur_linekey), SEEK_CUR);
                fwrite(&cur_linekey, sizeof(cur_linekey), 1, fl);
                /* Update the list in C */
                return_val = add_new_c_entry(info, cur_linekey.start_offset);
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
            fclose(fc);
            fc = NULL;
        }
        fwrite(&cur_linekey, sizeof(cur_linekey), 1, fl);
        //return_val = add_new_c_entry(info, cur_linekey.start_offset);
        goto out;
    }

  out:
    if (fl) fclose(fl);
    fl = NULL;
    if (fc) fclose(fc);
    fc = NULL;
    return return_val;
}

int search_info(linekey_info_t *info_array, const int lk_idx) {
    int return_val = 0;
    FILE *f = NULL;
    linekey_info_t cur_info;

    f = fopen(PATH_TO_C, "rb+");
    if (!f) {
        return -1;
    }

    /* No heder for now, could change in the future */
    /* Find the list for the linekey at position lk_idx */
    fseek(f, sizeof(linekey_info_t)*lk_idx, SEEK_SET);
    

  out:
    return return_val;
}

void print_linekey(linekey_t linekey) {
    printf("Value: 0x%016" PRIx64 " :: Offset: %lu :: Count: %lu\n", linekey.value, linekey.start_offset, linekey.count);
}

void print_info(linekey_info_t info, int cur_loc) {
    printf("ID: %d :: Position: %d :: Next: %d :: Location in file: %d\n", info.id, info.p, info.next_offset, cur_loc);
}

int print_data(const char *name) {
    int return_val = 0;
    FILE *f = NULL;
    int type = 0; /* 0 -> linekey || 1 -> info */
    linekey_t cur_linekey;
    linekey_info_t cur_info;

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
    } else {
        return_val = -1;
        goto out;
    }

    /* Read File -> currently print */
    if (type) {
        while (fread(&cur_info, sizeof(cur_info), 1, f) == 1)
            print_info(cur_info, ftell(f) - sizeof(cur_info));
    } else {
        while (fread(&cur_linekey, sizeof(cur_linekey), 1, f) == 1)
            print_linekey(cur_linekey);
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
    return 0;
}