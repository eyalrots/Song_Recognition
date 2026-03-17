#ifndef __DATABASE_H__
#define __DATABASE_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>

#define PATH_TO_L "../database/linekeys.bin"
#define L_NAME "LKDB"
#define PATH_TO_C "../database/song_info.bin"

typedef struct linekey {
    uint64_t value; // The value of the linekey
    uint64_t start_offset; // Where its list starts in C file
    uint64_t count; //Numebr of appearences of this linekey overall
} linekey_t;

typedef struct linekey_info {
    uint32_t id; // Song id in which linekey was found
    uint32_t p;  // Position in said song the linekey if from
    uint32_t next_offset; // Offset to the next item on the list
} linekey_info_t;

#endif