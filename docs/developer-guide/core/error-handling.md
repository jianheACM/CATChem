# Error Handling Developer Guide

CATChem is equipped with a comprehensive error handling system, defined in `src/core/error_mod.F90`, that provides robust error management, graceful degradation, and detailed error reporting. This guide describes the error handling system and how to use it effectively when developing for CATChem.

## Overview

The error handling system in `error_mod.F90` is a hybrid system that includes both a modern, object-oriented error management framework and a set of simpler, legacy error reporting functions.

- **Modern System (`ErrorManagerType`)**: A powerful, structured system for handling errors with context tracking, severity levels, and error categories. This is the recommended approach for new code.
- **Legacy System (`CC_Error`, `CC_Warning`)**: A set of simple subroutines for reporting errors and warnings. These are maintained for backward compatibility and are used by the modern system internally.

## The Modern Error Handling System: `ErrorManagerType`

The `ErrorManagerType` is the core of the modern error handling system. It provides a structured and powerful way to manage and report errors.

### Features

- **Error Context Stack**: The `ErrorManagerType` can maintain a stack of "contexts" (`ErrorContextType`). This allows for detailed error reports that show the call stack leading up to an error, which is invaluable for debugging.
- **Structured Error Information**: Errors are represented by an `ErrorInfoType` object, which stores a code, message, severity, category, and location.
- **Standardized Codes, Severities, and Categories**: The module provides a rich set of predefined error codes (e.g., `ERROR_INVALID_INPUT`, `ERROR_FILE_NOT_FOUND`), severity levels (e.g., `SEVERITY_WARNING`, `SEVERITY_ERROR`, `SEVERITY_FATAL`), and categories (e.g., `CATEGORY_INPUT`, `CATEGORY_COMPUTATION`).
- **Centralized Reporting**: The `report_error` procedure is the central point for reporting all errors.

### Using the `ErrorManagerType`

Here is a typical workflow for using the `ErrorManagerType` in a subroutine:

```fortran
subroutine my_data_processing_subroutine(data, error_mgr, rc)
  use error_mod
  implicit none

  real, intent(in) :: data(:)
  type(ErrorManagerType), intent(inout) :: error_mgr
  integer, intent(out) :: rc

  rc = CC_SUCCESS

  ! 1. Push the current context onto the stack
  call error_mgr%push_context("my_data_processing_subroutine", "Processing input data")

  ! 2. Perform operations and check for errors
  if (any(data < 0.0)) then
    ! 3. Report an error if one occurs
    call error_mgr%report_error(ERROR_INVALID_INPUT, "Input data cannot be negative", rc)
    call error_mgr%pop_context() ! Pop context before returning
    return
  end if

  ! ... process data ...

  ! 4. Pop the context from the stack before exiting
  call error_mgr%pop_context()

end subroutine my_data_processing_subroutine
```

#### 1. Pushing Context

Before performing any operations that might fail, you should push the current context onto the error manager's stack using `push_context`. This provides valuable information if an error occurs within the subroutine or in any of the subroutines it calls.

#### 2. Reporting Errors

If an error is detected, you should report it using the `report_error` procedure. You provide an error code, a descriptive message, and the return code variable.

#### 3. Popping Context

It is crucial to pop the context from the stack using `pop_context` before the subroutine exits. This should be done for all exit paths, including error returns.

### Error Codes, Severities, and Categories

The `error_mod.F90` module defines a wide range of error codes, severities, and categories. When reporting an error, you should choose the most appropriate code. The `ErrorManagerType` will automatically determine the severity and category from the code.

- **Error Codes**: e.g., `ERROR_INVALID_CONFIG`, `ERROR_MEMORY_ALLOCATION`, `ERROR_CONVERGENCE`.
- **Severity Levels**: `SEVERITY_INFO`, `SEVERITY_WARNING`, `SEVERITY_ERROR`, `SEVERITY_CRITICAL`, `SEVERITY_FATAL`. The error manager can be configured to abort execution on critical or fatal errors.
- **Error Categories**: `CATEGORY_INPUT`, `CATEGORY_COMPUTATION`, `CATEGORY_MEMORY`, `CATEGORY_IO`, etc. These are used for statistical tracking of errors.

## The Legacy Error Handling System

For simpler cases, or in code that has not yet been updated to use the `ErrorManagerType`, you can use the legacy `CC_Error` and `CC_Warning` subroutines.

### `CC_Error`

The `CC_Error` subroutine prints a formatted error message to standard output and sets the return code to `CC_FAILURE`.

```fortran
! Example of using CC_Error
subroutine legacy_error_example(temperature, rc)
  use error_mod, only: CC_Error, CC_SUCCESS, CC_FAILURE
  implicit none

  real, intent(in) :: temperature
  integer, intent(out) :: rc

  rc = CC_SUCCESS

  if (temperature < 0.0) then
    call CC_Error("Temperature cannot be negative", rc, "legacy_error_example")
    return
  end if

end subroutine legacy_error_example
```

### `CC_Warning`

The `CC_Warning` subroutine is similar to `CC_Error`, but it is for non-critical issues. It prints a warning message but does not change the return code.

```fortran
! Example of using CC_Warning
if (pressure > 110000.0) then
  call CC_Warning("Pressure is unusually high", rc, "pressure_check")
end if
```

## Best Practices for Error Handling

- **Use `ErrorManagerType` for New Code**: All new code should use the modern `ErrorManagerType` for structured and detailed error handling.
- **Push and Pop Context**: Always push the context at the beginning of a subroutine and pop it before every exit point.
- **Be Specific in Error Messages**: Provide clear, descriptive error messages that will help a user or developer diagnose the problem.
- **Choose the Right Error Code**: Use the most specific error code available from the list in `error_mod.F90`.
- **Check Return Codes**: Always check the return code of any subroutine that can fail.
- **Handle Errors Gracefully**: Whenever possible, your code should handle errors gracefully, clean up any allocated resources, and return an informative error code.
