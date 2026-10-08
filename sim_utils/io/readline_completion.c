#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <readline/readline.h>

//+
// Tab completion support for programs using read_a_line (see read_line.c).
//
// A program (e.g. Tao) registers a completion callback via readline_set_completion_fn.
// The callback is handed the current line buffer and cursor position and fills a buffer
// with newline-separated lines: the first line is the common prefix readline should
// insert, and each following line is one candidate.
//
// Callback return value:
//   n > 0   -- Number of candidates in buf (after the prefix line).
//   n == 0  -- Recognized completion position but nothing matches. Nothing is offered
//              and filename completion is suppressed.
//   n == -1 -- Not the program's to complete: fall back to readline's default
//              filename completion.
//
// Registration saves the previous readline completion settings so an embedder can
// restore them with readline_clear_completion_fn. If no callback is registered,
// readline behaves as before.
//
// Limitation: words break on blanks and tabs only (see readline_set_completion_fn),
// and that break set also governs the file name fallback. A file name glued to
// shell syntax, as in "spawn cat <./di<TAB>", is completed as the whole word "<./di"
// and finds nothing; put a blank before the name. Quoted names such as
// read "my dir/fi<TAB> do complete since the word starts after the quote.
//-

typedef int (*sim_rl_complete_fn)(const char* line, int point, int start, int end,
                                  char* buf, int buf_size);

static sim_rl_complete_fn completion_fn = NULL;

static int saved_state = 0;
static rl_completion_func_t* prev_completion_function = NULL;
static const char* prev_word_break_characters = NULL;
static const char* prev_completer_quote_characters = NULL;

#define SIM_RL_BUF_SIZE 65536
static char candidate_buf[SIM_RL_BUF_SIZE];

//----------------------------------------------------------------------------

static char** sim_rl_attempted_completion(const char* text, int start, int end) {
  (void) text;

  if (completion_fn == NULL) return NULL;

  int n = completion_fn(rl_line_buffer, rl_point, start, end, candidate_buf, SIM_RL_BUF_SIZE);

  if (n < 0) return NULL;   // Not ours: let readline do filename completion.

  rl_attempted_completion_over = 1;
  if (n == 0) return NULL;

  // matches[0] = common prefix (buffer line 1), matches[1..n] = candidates,
  // matches[n+1] = NULL. readline takes ownership and frees every string plus
  // the array itself, so everything is freshly malloc'd.

  char** matches = (char**) malloc((n + 2) * sizeof(char*));
  if (!matches) return NULL;

  char* p = candidate_buf;
  int i, n_filled = 0;   // Number of matches[] entries holding malloc'd strings.

  for (i = 0; i <= n && *p; i++) {
    char* nl = strchr(p, '\n');
    size_t len = nl ? (size_t)(nl - p) : strlen(p);

    char* str = (char*) malloc(len + 1);
    if (!str) break;
    memcpy(str, p, len);
    str[len] = '\0';
    matches[i] = str;
    n_filled++;

    if (!nl) break;
    p = nl + 1;
  }

  if (n_filled < 2) {   // Need at least the prefix and one candidate.
    for (i = 0; i < n_filled; i++) free(matches[i]);
    free(matches);
    return NULL;
  }

  matches[n_filled] = NULL;

  // A sole candidate ending in a prefix character (eg "quadrupole::" or "2@") is
  // meant to be typed onward, so do not append the usual trailing space.
  if (n_filled == 2) {
    size_t len = strlen(matches[1]);
    if (len > 0 && (matches[1][len-1] == ':' || matches[1][len-1] == '@')) {
      rl_completion_suppress_append = 1;
    }
  }

  return matches;
}

//----------------------------------------------------------------------------
//+
// Routine readline_set_completion_fn (sim_rl_complete_fn fn)
//
// Register a completion callback and install the completion hooks into readline,
// saving the previous settings on first use. Called from Fortran via a bind(c)
// interface (see tao_completion_mod.f90).
//-

void readline_set_completion_fn(sim_rl_complete_fn fn) {
  if (!saved_state) {
    prev_completion_function = rl_attempted_completion_function;
    prev_word_break_characters = rl_completer_word_break_characters;
    prev_completer_quote_characters = rl_completer_quote_characters;
    saved_state = 1;
  }
  completion_fn = fn;
  rl_attempted_completion_function = sim_rl_attempted_completion;
  // Tokens break on whitespace only so constructs like "2@q1", "orbit.x", and
  // "-universe" complete as single words. Must agree with the Fortran engine.
  rl_completer_word_break_characters = " \t";
  // Inside an unclosed quote the word starts after the quote, so quoted file names
  // with blanks complete on the fallback path. Program tokens are never quoted: the
  // engine sees the word start differ from its own and offers nothing, as before.
  rl_completer_quote_characters = "\"'";
}

//----------------------------------------------------------------------------
//+
// Routine readline_clear_completion_fn ()
//
// Remove the completion callback and restore the readline settings that were in
// effect before readline_set_completion_fn was first called.
//-

void readline_clear_completion_fn(void) {
  completion_fn = NULL;
  if (saved_state) {
    rl_attempted_completion_function = prev_completion_function;
    // Cast: the header declares this char* in some readline versions, const char* in others.
    rl_completer_word_break_characters = (char*) prev_word_break_characters;
    rl_completer_quote_characters = (char*) prev_completer_quote_characters;
  }
}
