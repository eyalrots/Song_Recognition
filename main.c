#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "func.h"

#define PATH_TO_FILE "../../sampled_audio/output.wav"

int main()
{
	unsigned char *data_buffer = NULL;
	double *audio_data = NULL;
	linekey_t *unique_linekeys = NULL;
	int unique_len = 0;
	long size = 0;
	wav_header_t header;
	int song_idx = 0;
	int i = 0;
	uint64_t cur_linekey = 0;

	int read_num = read_file(PATH_TO_FILE, &data_buffer, &size, &header);
	if (!data_buffer || read_num != 0) {
		fprintf(stderr, "Error getting data from file.\n");
		return 1;
	}

	convert_to_samples(data_buffer, &audio_data, &size,
			   header.bits_per_sample);

	// prepare arrays
	unique_linekeys = (linekey_t *)calloc(500, sizeof(linekey_t));
	unique_len = 500;
	for (i = 0; i < unique_len; i++) {
		unique_linekeys[i].position_arr =
			(int **)calloc(20, sizeof(int *));
		unique_linekeys[i].pos_arr_len = 20;
		unique_linekeys[i].new_pos = 0;
	}

	analyze_data(audio_data, size, header, song_idx, unique_linekeys,
		     &unique_len);

	printf("Final length of unique array is: %d\n", unique_len);

	for (i = 0; i < unique_len; i++) {
		cur_linekey = convert_linekey_to_number(
			unique_linekeys[i].binary_values);
		printf("linekey[%d] = %llu\n", i, cur_linekey);
	}

	printf("Done generation of linekeys, please continue to making the database.\n");

	free(unique_linekeys);
	free(data_buffer);
	free(audio_data);
	return 0;
}