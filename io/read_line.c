#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <readline/readline.h>
#include <readline/history.h>

//+
// Routine read_line_(char* tag, char* str, char* hist_file, int tag_len, int str_len, int hist_len)
//
// Routine to interface between the readline routine and the Fortran routine read_a_line.
// See read_a_line for more details.
//-

void read_line_(char* tag, char* str, char* hist_file, int tag_len, int str_len, int hist_len) {
  /* printf ("\n"); */
  char* str2 = readline (tag);

  if (!str2) {      // cntl-D pressed
    int ii;
    for (ii = 1; ii < str_len; ii++) 
      str[ii] = ' ';
    str[0] = 24;
    return;
  }

  strcpy (str, str2);
  int ii;
  for (ii = strlen(str2); ii < str_len; ii++) 
    str[ii] = ' ';

  /* If the line has any text in it, save it on the history stack */

  if (str2 && *str2) {
    add_history (str2);

    if (hist_len > 1) {
      if (hist_file[0] == '~') {
        char* home = getenv("HOME");
        const size_t home_len = strlen(home);
        char* name = malloc(hist_len + home_len - 1);
        memcpy (name, home, home_len);
        memcpy (name + home_len, hist_file+1, hist_len-1);  // copy with out leading "~".
        int stat = append_history(1, name);
        free(name);

      } else {
        int stat = append_history(1, hist_file);
      }
    }
  }

  free (str2);
}

//+
// Routine read_history_(char* file, int* stat, int file_len)
//
// Routine to interface between the readline routine read_history and the Fortran routine readline_read_history.
//-

void read_history_(char* file, int* stat, int file_len) {
  if (file[0] == '~') {
    char* home = getenv("HOME");
    const size_t len_h = strlen(home);
    const size_t len_f = strlen(file);
    char* name = malloc(len_f + len_h);
    memcpy (name, home, len_h);
    memcpy (name + len_h, file+1, len_f);  // copy with out leading "~".
    *stat = read_history(name);
    free(name);

  } else {
    *stat = read_history(file);
  }
}

//+
// Routine write_history_(char* file, int* stat, int file_len)
//
// Routine to interface between the readline routine write_history and the Fortran routine readline_write_history.
//-

void write_history_(char* file, int* stat, int file_len) {
  if (file[0] == '~') {
    char* home = getenv("HOME");
    const size_t len_h = strlen(home);
    const size_t len_f = strlen(file);
    char* name = malloc(len_f + len_h);
    memcpy (name, home, len_h);
    memcpy (name + len_h, file+1, len_f);  // copy with out leading "~".
    *stat = write_history(name);
    free(name);

  } else {
    *stat = write_history(file);
  }
}
