
#' @include classes.R
NULL

.name = NULL # dummy

#' @name g_assert
#' @title Make guard clause assertions
#' @param x object to make assertions of
#' @param test assertion test to perform on \code{x}
#' @param n Frames back in which to evaluate the gobject param
#' @param msg error message provided as character vector. Using the token .name
#' will send the name of object being tested in the stack frame referenced by
#' \code{n}
#' @param ... additional params to pass
#' @keywords internal
#' @noRd
g_assert = function(x, test, msg = NULL, n = 2L, ...) {
  if(!test) {
    # get name of function where test failed
    fn_name = deparse(sys.calls()[[sys.nframe() - n]])
    # get name of object that failed test
    .name = deparse(eval(call('substitute', as.name(substitute(x)), parent.frame(n = 1L))))
    .name = paste0('\"\'', .name, '\'\"')

    # compose message
    msg = gsub(pattern = '.name', replacement = .name, x = deparse(substitute(msg)))
    msg = parse(text = msg)

    # send error
    stop(wrap_txt(fn_name, ':\n', eval(msg, envir = parent.frame(n = 2L)), errWidth = TRUE),
         call. = FALSE)
  }
}




#' @param gobject giotto object
#' @param n Frames back in which to evaluate the gobject param
#' @keywords internal
#' @noRd
guard_against_notgiotto = function(gobject, n = 1L, ...) {
  fn_name = deparse(sys.calls()[[sys.nframe() - n]])
  orig_name = deparse(eval(call('substitute', as.name(substitute(gobject)), parent.frame())))
  if(!hasArg(gobject)) stop(wrap_txt(fn_name, ':\ngiotto object must be given',
                                     errWidth = TRUE),
                            call. = FALSE)
  if(!inherits(gobject, 'giotto')) {
    stop(wrap_txt(fn_name, ':\n', orig_name, 'is not a giotto object',
                  errWidth = TRUE),
         call. = FALSE)
  }
}




assert_DT = function(x) {
  g_assert(
    x, test = inherits(x, 'data.table'),
    msg = c(.name, 'must be of class data.table, not', class(x))
  )
}



assert_file = function(x) {
  g_assert(
    x, test = is.character(x),
    msg = c(.name, 'must be a character vector filepath')
  )
  g_assert(
    x, test = file.exists(x),
    msg = c(.name, 'is not an existing file')
  )
}



assert_numeric = function(x) {
  g_assert(
    x, test = is.numeric(x),
    msg = c(.name, 'must be of class numeric, not', class(x))
  )
}


