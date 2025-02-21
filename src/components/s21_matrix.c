#include "../s21_matrix.h"

int is_matrix_ok(matrix_t *A) {
  if (A != NULL && A->matrix != NULL && A->rows >= 1 && A->columns >= 1)
    return SUCCESS;
  else
    return FAILURE;
}

int same_size(matrix_t *A, matrix_t *B) {
  if (A->rows == B->rows && A->columns == B->columns)
    return SUCCESS;
  else
    return FAILURE;
}

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  if (result == NULL || rows < 1 || columns < 1) return ERR_INCORRECT_MATRIX;

  int status = OK;
  result->rows = rows;
  result->columns = columns;
  result->matrix = (double **)calloc(rows, sizeof(double *));
  if (result->matrix == NULL) {
    status = ERR_INCORRECT_MATRIX;
  } else {
    int i = 0;
    while (i < rows && status == OK) {
      result->matrix[i] = (double *)calloc(columns, sizeof(double));
      if (result->matrix[i] == NULL) {
        for (int j = 0; j < i; j++) {
          free(result->matrix[j]);
        }
        free(result->matrix);
        status = ERR_INCORRECT_MATRIX;
      }
      i++;
    }
  }
  return status;
}

void s21_remove_matrix(matrix_t *A) {
  if (is_matrix_ok(A)) {
    for (int i = 0; i < A->rows; i++) {
      free(A->matrix[i]);
      A->matrix[i] = NULL;
    }
    free(A->matrix);
    A->matrix = NULL;
  }
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  if (!is_matrix_ok(A) || !is_matrix_ok(B) || !same_size(A, B)) return FAILURE;

  int status = SUCCESS;
  for (int i = 0; status == SUCCESS && i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      if (fabs(A->matrix[i][j] - B->matrix[i][j]) >= 1e-7) status = FAILURE;
    }
  }
  return status;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  if (!is_matrix_ok(A) || !is_matrix_ok(B) || result == NULL)
    return ERR_INCORRECT_MATRIX;
  if (!same_size(A, B)) return CALCULATION_ERROR;

  int status = OK;

  if ((status = s21_create_matrix(A->rows, A->columns, result)) == OK) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++)
        result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
    }
  }
  return status;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  if (!is_matrix_ok(A) || !is_matrix_ok(B) || result == NULL)
    return ERR_INCORRECT_MATRIX;
  if (!same_size(A, B)) return CALCULATION_ERROR;

  int status = OK;
  if ((status = s21_create_matrix(A->rows, A->columns, result)) == OK) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++)
        result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
    }
  }
  return status;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  if (!is_matrix_ok(A) || result == NULL) return ERR_INCORRECT_MATRIX;

  int status = OK;
  if ((status = s21_create_matrix(A->rows, A->columns, result)) == OK) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++)
        result->matrix[i][j] = A->matrix[i][j] * number;
    }
  }

  return status;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  if (!is_matrix_ok(A) || !is_matrix_ok(B) || (result == NULL))
    return ERR_INCORRECT_MATRIX;
  if (A->columns != B->rows) return CALCULATION_ERROR;

  int status = OK;

  if ((status = s21_create_matrix(A->rows, B->columns, result)) == OK) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < B->columns; j++) {
        for (int k = 0; k < A->columns; k++) {
          result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
        }
      }
    }
  }
  return status;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  if (!is_matrix_ok(A) || result == NULL) return ERR_INCORRECT_MATRIX;
  int status = OK;
  if ((status = s21_create_matrix(A->columns, A->rows, result)) == OK) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++)
        result->matrix[j][i] = A->matrix[i][j];
    }
  }
  return status;
}

int create_minor_matrix(matrix_t *A, int exclude_row, int exclude_col,
                        matrix_t *submatrix) {
  if (!is_matrix_ok(A) || submatrix == NULL) return ERR_INCORRECT_MATRIX;
  int status = OK;
  if ((status = s21_create_matrix(A->rows - 1, A->columns - 1, submatrix)) ==
      OK) {
    int sub_i = 0;
    for (int i = 0; i < A->rows; i++) {
      int sub_j = 0;
      for (int j = 0; j < A->columns; j++) {
        int row_condition = (i != exclude_row);
        int col_condition = (j != exclude_col);
        if (row_condition && col_condition) {
          submatrix->matrix[sub_i][sub_j] = A->matrix[i][j];
          sub_j++;
        }
      }
      if (i != exclude_row) {
        sub_i++;
      }
    }
  }
  return status;
}

int s21_determinant(matrix_t *A, double *result) {
  if (!is_matrix_ok(A) || result == NULL) return ERR_INCORRECT_MATRIX;
  if (A->rows != A->columns) return CALCULATION_ERROR;

  int status = OK;
  if (A->rows == 1)
    *result = A->matrix[0][0];
  else if (A->rows == 2)
    *result =
        A->matrix[0][0] * A->matrix[1][1] - A->matrix[0][1] * A->matrix[1][0];
  else {
    int exit_flag = 0;
    int sign = 1;
    *result = 0;
    for (int i = 0; i < A->columns && !exit_flag; i++) {
      matrix_t submatrix = {0};
      if ((status = create_minor_matrix(A, 0, i, &submatrix)) != OK)
        exit_flag = 1;

      double determinant = 0.0;
      if ((status = s21_determinant(&submatrix, &determinant)) != OK) {
        s21_remove_matrix(&submatrix);
        exit_flag = 1;
      }
      *result += sign * A->matrix[0][i] * determinant;
      sign = -sign;
      s21_remove_matrix(&submatrix);
    }
  }
  return status;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  if (!is_matrix_ok(A) || result == NULL) return ERR_INCORRECT_MATRIX;
  if (A->rows != A->columns) return CALCULATION_ERROR;
  int status = OK;
  if ((status = s21_create_matrix(A->rows, A->columns, result)) == OK) {
    if (result->rows == 1 && result->columns == 1)
      result->matrix[0][0] = A->matrix[0][0];
    else {
      int exit_flag = 0;
      for (int i = 0; i < A->rows && !exit_flag; i++) {
        for (int j = 0; j < A->columns; j++) {
          matrix_t minor;
          double determinant = 0.0;
          exit_flag = ((status = create_minor_matrix(A, i, j, &minor))) != OK
                          ? 1
                          : exit_flag;
          exit_flag = ((status = s21_determinant(&minor, &determinant))) != OK
              ? (s21_remove_matrix(&minor)),
          1 : exit_flag;
          result->matrix[i][j] = determinant * ((i + j) % 2 == 0 ? 1 : -1);
          s21_remove_matrix(&minor);
        }
      }
    }
  }
  return status;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  if (!is_matrix_ok(A) || result == NULL) return ERR_INCORRECT_MATRIX;

  if (A->rows != A->columns) return CALCULATION_ERROR;

  int status = OK;
  double determinant = 0.0;

  if ((status = s21_determinant(A, &determinant)) == OK) {
    if (fabs(determinant) <= 1e-7) status = CALCULATION_ERROR;
  }

  if (status == OK) {
    if (A->rows == 1) {
      if ((status = s21_create_matrix(A->rows, A->columns, result)) == OK)
        result->matrix[0][0] = 1 / A->matrix[0][0];
    } else {
      matrix_t temp = {0};
      matrix_t transposed = {0};

      status = s21_calc_complements(A, &temp);

      if (status == OK) status = s21_transpose(&temp, &transposed);

      if (status == OK)
        status = s21_mult_number(&transposed, 1 / determinant, result);

      s21_remove_matrix(&temp);
      s21_remove_matrix(&transposed);
    }
  }
  return status;
}