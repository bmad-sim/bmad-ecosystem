#include <stdlib.h>
#include <string.h>
#include <readline/readline.h>

//+
// Tab completion support for programs using read_a_line (see read_line.c).
//
// A program (e.g. Tao) registers a completion callback via readline_set_completion_fn.
// The callback is handed the current line buffer and cursor position and fills a buffer
// with newline-separated completion candidates.
//
// Callback return value:
//   n >= 0  -- Number of candidates placed in buf. Filename completion is suppressed
//              (even when n == 0).
//   n == -1 -- The word being completed is a file path: fall back to readline's
//              default filename completion.
//
// If no callback is registered, readline behaves as before (filename completion).
//-

typedef int (*sim_rl_complete_fn)(const char* line, int point, int start, int end,
                                  char* buf, int buf_size);

static sim_rl_complete_fn completion_fn = NULL;

#define SIM_RL_BUF_SIZE 65536
static char candidate_buf[SIM_RL_BUF_SIZE];

//----------------------------------------------------------------------------

static char** sim_rl_attempted_completion(const char* text, int start, int end) {
  (void) text;

  if (completion_fn == NULL) return NULL;

  int n = completion_fn(rl_line_buffer, rl_point, start, end, candidate_buf, SIM_RL_BUF_SIZE);

  if (n < 0) return NULL;   // File path: let readline do filename completion.

  rl_attempted_completion_over = 1;   // Never fall back to filename completion.
  if (n == 0) return NULL;

  // matches[0] = longest common prefix, matches[1..n] = candidates, matches[n+1] = NULL.
  // readline takes ownership and frees every string plus the array itself.

  char** matches = (char**) malloc((n + 2) * sizeof(char*));
  if (!matches) return NULL;

  char* p = candidate_buf;
  int i, n_found = 0;
  size_t lcp_len = 0;

  for (i = 1; i <= n && *p; i++) {
    char* nl = strchr(p, '\n');
    size_t len = nl ? (size_t)(nl - p) : strlen(p);

    char* cand = (char*) malloc(len + 1);
    if (!cand) break;
    memcpy(cand, p, len);
    cand[len] = '\0';
    matches[n_found + 1] = cand;

    if (n_found == 0) {
      lcp_len = len;
    } else {
      size_t j = 0;
      while (j < lcp_len && j < len && matches[1][j] == cand[j]) j++;
      lcp_len = j;
    }
    n_found++;

    if (!nl) break;
    p = nl + 1;
  }

  if (n_found == 0) {
    free(matches);
    return NULL;
  }

  matches[0] = (char*) malloc(lcp_len + 1);
  if (!matches[0]) {
    for (i = 1; i <= n_found; i++) free(matches[i]);
    free(matches);
    return NULL;
  }
  memcpy(matches[0], matches[1], lcp_len);
  matches[0][lcp_len] = '\0';

  matches[n_found + 1] = NULL;
  return matches;
}

//----------------------------------------------------------------------------
//+
// Routine readline_set_completion_fn (sim_rl_complete_fn fn)
//
// Register a completion callback and install the completion hooks into readline.
// Called from Fortran via a bind(c) interface (see tao_completion_mod.f90).
//-

void readline_set_completion_fn(sim_rl_complete_fn fn) {
  completion_fn = fn;
  rl_attempted_completion_function = sim_rl_attempted_completion;
  // Tokens break on whitespace only so constructs like "2@q1", "orbit.x", and
  // "-universe" complete as single words. Must agree with the Fortran engine.
  rl_completer_word_break_characters = " \t";
}
