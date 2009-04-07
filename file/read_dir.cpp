#include "CESR_platform.h"

#include <dirent.h>
#include <stdlib.h>
#include <string.h>

DIR *ptr_dir;
struct dirent* dir_ent = NULL;


//-----------------------------------------------------------------------

#if defined(CESR_VMS)

extern "C" void set_string(const char*, const int&, char*, const int&);

extern "C" void open_dir(char* dir, bool& valid) {
  if (ptr_dir != NULL) closedir(ptr_dir);
  ptr_dir = opendir(dir);
  valid = (ptr_dir != NULL);
}

extern "C" void close_dir() {
  if (ptr_dir) closedir(ptr_dir);
}

extern "C" void read_dir(char* file, const int& n_len, bool& valid) {
  dir_ent = readdir(ptr_dir);
  if (dir_ent == NULL) {
    closedir(ptr_dir);
    valid = false;
    return;
  }

  int str_len = strlen(dir_ent->d_name);
  set_string(dir_ent->d_name, str_len, file, n_len);
  valid = true;
  return;

}

//-----------------------------------------------------------------------

#elif defined(CESR_WINCVF)

extern "C" void set_string_(const char*, const int&, char*, const int&);

extern "C" void CLOSE_DIR() {
  if (ptr_dir) closedir(ptr_dir);
  ptr_dir = NULL;
}

extern "C" void OPEN_DIR(char* dir, bool& valid) {
  if (ptr_dir != NULL) close_dir_();
  ptr_dir = opendir(dir);
  valid = (ptr_dir != NULL);
}

extern "C" void READ_DIR(char* file, const int& n_len, bool& valid) {
  dir_ent = readdir(ptr_dir);
  if (dir_ent == NULL) {
    close_dir_();
    valid = false;
    return;
  }

  int str_len = strlen(dir_ent->d_name);
  set_string_(dir_ent->d_name, str_len, file, n_len);
  valid = true;
  return;

}

//-----------------------------------------------------------------------

#else

extern "C" void set_string_(const char*, const int&, char*, const int&);

extern "C" void close_dir_() {
  if (ptr_dir) closedir(ptr_dir);
  ptr_dir = NULL;
}

extern "C" void open_dir_(char* dir, bool& valid) {
  if (ptr_dir != NULL) close_dir_();
  ptr_dir = opendir(dir);
  valid = (ptr_dir != NULL);
}

extern "C" void read_dir_(char* file, const int& n_len, bool& valid) {
  dir_ent = readdir(ptr_dir);
  if (dir_ent == NULL) {
    close_dir_();
    valid = false;
    return;
  }

  int str_len = strlen(dir_ent->d_name);
  set_string_(dir_ent->d_name, str_len, file, n_len);
  valid = true;
  return;

}



#endif
