BEGIN{
	IGNORECASE = 1;

  n_modules = 0;
}

/module/ {
	if (tolower($1) == "module") {
    if (tolower($2) != "procedure") {
      module_name[n_modules++] = tolower($2);
    }
  }
}

END {
  n_fields = split(FILENAME, name_fields, /\//);
  filename = name_fields[n_fields];
#  printf("Filename is:  %s\n", filename);

  for (i = 0; i < n_modules; i++) {
    printf("%s/%s.mod: %s\n", OUT_DIR, module_name[i], filename);
    printf("\t%s\n",RULE);
  }
}
