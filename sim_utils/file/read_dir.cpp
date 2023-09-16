#include <dirent.h>
#include <stdlib.h>
#include <string.h>

DIR *ptr_dir;
struct dirent* dir_ent = NULL;


extern "C" void set_fortran_string_(const char*, const int&, char*, const int&);

extern "C" void close_dir_() {
  if (ptr_dir) closedir(ptr_dir);
  ptr_dir = NULL;
}

extern "C" void open_dir_(char* dir, int& valid) {
  if (ptr_dir != NULL) close_dir_();
  ptr_dir = opendir(dir);
  valid = (ptr_dir != NULL);
}

extern "C" void read_dir_(char* file, const int& n_len, int& valid) {
  dir_ent = readdir(ptr_dir);
  if (dir_ent == NULL) {
    close_dir_();
    valid = false;
    return;
  }

  int str_len = strlen(dir_ent->d_name);
  set_fortran_string_(dir_ent->d_name, str_len, file, n_len);
  valid = true;
  return;

}

