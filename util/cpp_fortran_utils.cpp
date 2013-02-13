extern "C" void bool_to_int (bool& bool_logic, int& int_logic) {
  int_logic = bool_logic;
}

extern "C" void int_to_bool (int& int_logic, bool& bool_logic) {
  if (int_logic == 0) {
    bool_logic = false;
  } else {
    bool_logic = true;
  }
}
