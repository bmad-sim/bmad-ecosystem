BEGIN{
	IGNORECASE = 1;
}

/module/ {
	if (tolower($1) == "module") {
    if (tolower($2) != "procedure") {
      printf("%s ", tolower($2));	
    }
  }
}
