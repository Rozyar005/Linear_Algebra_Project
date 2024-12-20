import java.util.Scanner;


public class Main {


        public static void main(String[] args) {

                Scanner scanner = new Scanner(System.in);
                System.out.println("Welcome to the Matrix Algebra Program!");

                int choice = 0;

                do {
                    System.out.println("Choose an operation:");
                    System.out.println("1. Matrix Addition");
                    System.out.println("2. Matrix Subtraction");
                    System.out.println("3. Matrix Multiplication");
                    System.out.println("4. Scalar Multiplication");
                    System.out.println("5. Transpose of a Matrix");
                    System.out.println("6. Determinant of a Matrix");
                    System.out.println("7. Inverse of a Matrix");
                    System.out.println("8. Solve Linear System of Equations");
                    System.out.println("9. Find Eigenvalues and Eigenvectors");
                    System.out.println("10. Gauss-Jordan or Cramer");
                    System.out.println("-1. Exit");
                    System.out.print("Enter your choice: ");
                    choice=scanner.nextInt();

                    // Process choice
                    switch (choice) {
                        case 1 -> matrixAddition(scanner);
                        case 2 -> matrixSubtraction(scanner);
                        case 3 -> matrixMultiplication(scanner);
                        case 4 -> scalarMultiplication(scanner);
                        case 5 -> transposeMatrix(scanner);
                        case 6 -> determinantMatrix(scanner);
                        case 7 -> inverseMatrix(scanner);
                        case 8 -> solveLinearSystem(scanner);
                        case 9 -> findEigenValuesAndVectors(scanner);
                        case 10 -> cal();
                        case -1 -> System.out.println("Goodbye!");
                        default -> System.out.println("Invalid choice! Please try again.");
                    }

                } while (choice != -1);
         }



        //This method is to input the values of the matrix each column's row
        private static double[][] inputMatrix(Scanner scanner, String name) {
            System.out.println("\nEnter dimensions of " + name + " matrix (rows and columns):");
            int rows = scanner.nextInt();
            int cols = scanner.nextInt();
            double[][] matrix = new double[rows][cols];
            System.out.println("Enter elements of " + name + " matrix:");
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    System.out.print("[" + (i+1) + "][" + (j+1) + "]: ");
                    matrix[i][j] = scanner.nextDouble();
                }
            }
            printMatrix(matrix, "You entered the following matrix:");
            return matrix;
        }

        //This is for printing the matrix
        private static void printMatrix(double[][] matrix, String message) {
            System.out.println(message);
            for (double[] row : matrix) {
                for (double val : row) {
                    System.out.printf("%.2f ", val);
                }
                System.out.println();
            }
        }
        //This is for the additionn
        private static void matrixAddition(Scanner scanner) {
            double[][] A = inputMatrix(scanner, "first");
            double[][] B = inputMatrix(scanner, "second");

            if (A.length != B.length || A[0].length != B[0].length) {
                System.out.println("Matrix dimensions do not match for addition.");
                return;
            }

            double[][] result = new double[A.length][A[0].length];
            for (int i = 0; i < A.length; i++) {    //In here it adds the matrixes first check teh matrix length then the matrix column
                for (int j = 0; j < A[0].length; j++) {  //This is for the matrix column
                    result[i][j] = A[i][j] + B[i][j];
                }
            }

            printMatrix(result, "Resultant Matrix after Addition:");
        }

        private static void matrixSubtraction(Scanner scanner) {
            double[][] A = inputMatrix(scanner, "first");
            double[][] B = inputMatrix(scanner, "second");

            if (A.length != B.length || A[0].length != B[0].length) {
                System.out.println("Matrix dimensions do not match for subtraction.");
                return;
            }

            double[][] result = new double[A.length][A[0].length];
            for (int i = 0; i < A.length; i++) { //same as the addition checks the row's length then the cloumn but for subtraction this time
                for (int j = 0; j < A[0].length; j++) {
                    result[i][j] = A[i][j] - B[i][j];
                }
            }

            printMatrix(result, "Resultant Matrix after Subtraction:");
        }

        private static void matrixMultiplication(Scanner scanner) {
            double[][] A = inputMatrix(scanner, "first");
            double[][] B = inputMatrix(scanner, "second");

            if (A[0].length != B.length) { //Check point
                System.out.println("Matrix dimensions do not allow multiplication.");
                return;
            }

            double[][] result = new double[A.length][B[0].length]; //creates an array of 2 dimentional array that goes by the no.row=no of row of matrix A and no.col=no of col of matrix B
            for (int i = 0; i < A.length; i++) { //go by matrix A  row length
                for (int j = 0; j < B[0].length; j++) { //go by matrix B column length
                    for (int k = 0; k < A[0].length; k++) { //go by matrix A column length
                        result[i][j] += A[i][k] * B[k][j];
                    }
                }
            }

            printMatrix(result, "Resultant Matrix after Multiplication:");
        }

        private static void scalarMultiplication(Scanner scanner) {
            double[][] A = inputMatrix(scanner, "input");
            System.out.print("Enter scalar value: ");
            double scalar = scanner.nextDouble();

            double[][] result = new double[A.length][A[0].length]; //It uses to scalar one matrix
            for (int i = 0; i < A.length; i++) { //goes by the matrix row length
                for (int j = 0; j < A[0].length; j++) { //go by matrix column length
                    result[i][j] = A[i][j] * scalar;
                }
            }

            printMatrix(result, "Resultant Matrix after Scalar Multiplication:");
        }

        private static void transposeMatrix(Scanner scanner) {
            double[][] A = inputMatrix(scanner, "input");
            double[][] result = new double[A[0].length][A.length]; //creates an array to input the value of the resultant matrix by matrix A column length for the row of the resultant matrix and
            //matrix A row length for the resultant matrix Column length

            for (int i = 0; i < A.length; i++) {
                for (int j = 0; j < A[0].length; j++) {
                    result[j][i] = A[i][j];
                }
            }

            printMatrix(result, "Transpose of the Matrix:");
        }

        private static void determinantMatrix(Scanner scanner) {
            double[][] A1 = inputMatrix(scanner, "square");
            if (A1.length != A1[0].length) { //Check point
                System.out.println("Matrix must be square to calculate determinant.");
                return;
            }

            double determinant = calculateDeterminant(A1);
            System.out.println("Determinant of the Matrix: " + determinant);
        }


        //This is a helper method to calculate determinant
        private static double calculateDeterminant(double[][] matrix) {
            int n = matrix.length;
            if (n == 1) return matrix[0][0];
            if (n == 2) return (matrix[0][0] * matrix[1][1]) - (matrix[0][1] * matrix[1][0]);

            double det = 0;
            for (int col = 0; col < n; col++) {
                det += Math.pow(-1, col) * matrix[0][col] * calculateDeterminant(minorMatrix(matrix, 0, col));
            }
            return det;
        }

        //This is a helper method to calculate minorMatrix
        private static double[][] minorMatrix(double[][] matrix, int row, int col) {
            int n = matrix.length;
            double[][] minor = new double[n - 1][n - 1];
            for (int i = 0, mi = 0; i < n; i++) {
                if (i == row) continue;
                for (int j = 0, mj = 0; j < n; j++) {
                    if (j == col) continue;
                    minor[mi][mj++] = matrix[i][j];
                }
                mi++;
            }
            return minor;
        }

        //This is a helper method to calculate inverse of a matrix
        private static void inverseMatrix(Scanner scanner) {
            double[][] matrix = inputMatrix(scanner, "square");
            if (matrix.length != matrix[0].length) {
                System.out.println("Matrix must be square to find its inverse.");
                return;
            }

            double determinant = calculateDeterminant(matrix);
            if (determinant == 0) {
                System.out.println("Matrix is singular and cannot be inverted (determinant = 0).\n");
                return;
            }

            double[][] adjugate = new double[matrix.length][matrix.length];
            for (int i = 0; i < matrix.length; i++) {
                for (int j = 0; j < matrix.length; j++) {
                    adjugate[j][i] = Math.pow(-1, i + j) * calculateDeterminant(minorMatrix(matrix, i, j));
                }
            }

            double[][] inverse = new double[matrix.length][matrix.length];
            for (int i = 0; i < matrix.length; i++) {
                for (int j = 0; j < matrix.length; j++) {
                    inverse[i][j] = adjugate[i][j] / determinant;
                }
            }

            printMatrix(inverse, "Inverse of the Matrix:");
        }


        //this is for solving linear system
        private static void solveLinearSystem(Scanner scanner) {
            System.out.println("\nSolving Linear System of Equations (Ax = b):");
            double[][] A = inputMatrix(scanner, "coefficient");
            double[][] b = inputMatrix(scanner, "constant column");

            if (A.length != A[0].length || A.length != b.length || b[0].length != 1) {
                System.out.println("Invalid dimensions. A must be square, and b must be a column vector.");
                return;
            }

            double determinant = calculateDeterminant(A);
            if (determinant == 0) {
                System.out.println("System has no unique solution (determinant = 0).\n");
                return;
            }

            double[][] inverse = new double[A.length][A.length];
            for (int i = 0; i < A.length; i++) {
                for (int j = 0; j < A.length; j++) {
                    inverse[j][i] = Math.pow(-1, i + j) * calculateDeterminant(minorMatrix(A, i, j)) / determinant;
                }
            }

            double[][] x = matrixMultiplication(inverse, b);
            printMatrix(x, "Solution Vector x:");
        }


        //this is for eigen value and eigen vectors
        private static void findEigenValuesAndVectors(Scanner scanner) {
            System.out.println(" square matrix (2x2):");
            // Create matrix
            double[][] matrix = new double[2][2];

            // Input matrix elements
            System.out.println("Enter matrix elements:");
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    System.out.printf("Enter element for position [%d][%d]: ", (i+1), (j+1));
                    matrix[i][j] = scanner.nextDouble();
                }
            }
            // Calculate and print eigenvalues
            calculate2x2Eigenvalues(matrix);

        }
        private static void calculate2x2Eigenvalues(double[][] matrix) {
            //This is for using Quadratic Formula
            double a = matrix[0][0];
            double b = matrix[0][1];
            double c = matrix[1][0];
            double d = matrix[1][1];

            // Calculate trace and determinant
            double trace = a + d;
            double determinant = a * d - b * c;

            // Calculate eigenvalues using quadratic formula
            double discriminant = trace * trace - 4 * determinant;

            if (discriminant >= 0) {
                // Calculate two eigenvalues
                double eigenvalue1 = (trace + Math.sqrt(discriminant)) / 2;
                double eigenvalue2 = (trace - Math.sqrt(discriminant)) / 2;

                System.out.println("Eigenvalue 1: " + eigenvalue1);
                System.out.println("Eigenvalue 2: " + eigenvalue2);

                // Calculate eigenvectors for both eigenvalues
                System.out.println("\nEigenvector for Eigenvalue 1:");
                findEigenvector(matrix, eigenvalue1);

                System.out.println("\nEigenvector for Eigenvalue 2:");
                findEigenvector(matrix, eigenvalue2);
            } else {
                System.out.println("Complex eigenvalues (not supported in this simple version)");
            }
        }

        private static double[] findEigenvector(double[][] matrix, double eigenvalue) {
            System.out.printf("Calculating Eigenvector for Eigenvalue: %.4f\n", eigenvalue);

            double[][] modifiedMatrix = {
                    {matrix[0][0] - eigenvalue, matrix[0][1]},
                    {matrix[1][0], matrix[1][1] - eigenvalue}
            };

            double x = 1; // Easy for use Values
            double y = 0;

            // Use the first row of the modified matrix to compute the eigenvector
            if (Math.abs(modifiedMatrix[0][1]) > 1e-9) { // Avoid division by zero
                y = -modifiedMatrix[0][0] / modifiedMatrix[0][1];
            } else if (Math.abs(modifiedMatrix[1][1]) > 1e-9) {
                y = -modifiedMatrix[1][0] / modifiedMatrix[1][1];
            } else {
                System.out.println("Matrix has a singular row, unable to calculate eigenvector directly.");
            }

            System.out.printf("Eigenvector for Eigenvalue %.4f: x = %.4f, y = %.4f\n", eigenvalue, x, y);
            return new double[]{x, y};
        }

        public static double[][] matrixMultiplication(double[][] A, double[][] B) {
            double[][] result = new double[A.length][B[0].length];
            for (int i = 0; i < A.length; i++) {
                for (int j = 0; j < B[0].length; j++) {
                    for (int k = 0; k < A[0].length; k++) {
                        result[i][j] += A[i][k] * B[k][j];
                    }
                }
            }
            return result;
        }


    public static void cal() {
        Scanner scanner = new Scanner(System.in);

        System.out.println("Welcome to the Matrix Solver!");
        System.out.print("Enter the size of the square matrix (2 or 3): ");
        int size = scanner.nextInt();

        if (size != 2 && size != 3) {
            System.out.println("Invalid size. Please restart and choose 2 or 3.");
            return;
        }

        double[][] matrix = new double[size][size + 1];

        System.out.println("Enter the coefficients of the matrix row by row (including the augmented column):");
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size + 1; j++) {
                if (j < size) {
                    System.out.print("Enter coefficient for x" + (j + 1) + " at row " + (i + 1) + ": ");
                } else {
                    System.out.print("Enter constant for row " + (i + 1) + ": ");
                }
                matrix[i][j] = scanner.nextDouble();
            }
        }

        System.out.println("Choose the solving method: ");
        System.out.println("1. Gauss-Jordan Method");
        System.out.println("2. Cramer's Rule");
        int method = scanner.nextInt();

        if (method == 1) {
            printMatrixx(matrix, true);
            solveUsingGaussJordan(matrix, size);
        } else if (method == 2) {
            printMatrixx(matrix, true);
            solveUsingCramer(matrix, size);
        } else {
            System.out.println("Invalid choice. Please restart.");
        }

    }

    private static void printMatrixx(double[][] matrix, boolean includeConstants) {
        System.out.println("The entered matrix is:");
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                System.out.print(matrix[i][j] + "\t");
            }
            System.out.println();
        }
    }


    //This is the method for handling the calculation
    //we will explain the exact logic behind it in the day of presenting it
    private static void solveUsingGaussJordan(double[][] matrix, int size) {
        for (int i = 0; i < size; i++) {
            double pivot = matrix[i][i];
            for (int j = 0; j < size + 1; j++) {
                matrix[i][j] /= pivot;
            }

            for (int k = 0; k < size; k++) {
                if (k != i) {
                    double factor = matrix[k][i];
                    for (int j = 0; j < size + 1; j++) {
                        matrix[k][j] -= factor * matrix[i][j];
                    }
                }
            }
        }

        System.out.println("Solution using Gauss-Jordan Method:");
        for (int i = 0; i < size; i++) {
            System.out.println("x" + (i + 1) + " = " + matrix[i][size]);
        }
    }

    //this method is to calculate the cramer
    private static void solveUsingCramer(double[][] matrix, int size) {
        double[][] coefficients = new double[size][size];
        double[] constants = new double[size];

        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                coefficients[i][j] = matrix[i][j];
            }
            constants[i] = matrix[i][size];
        }

        double determinant = calculateDeterminantt(coefficients);
        if (determinant == 0) {
            System.out.println("The system has no unique solution.");
            return;
        }

        double[] solutions = new double[size];
        for (int i = 0; i < size; i++) {
            double[][] tempMatrix = copyMatrix(coefficients);
            for (int j = 0; j < size; j++) {
                tempMatrix[j][i] = constants[j];
            }
            solutions[i] = calculateDeterminantt(tempMatrix) / determinant;
        }
        //Final Solution is in here
        System.out.println("Solution using Cramer's Rule:");
        for (int i = 0; i < size; i++) {
            System.out.println("x" + (i + 1) + " = " + solutions[i]);
        }


    }
    //this method is to calculate the determinant
    private static double calculateDeterminantt(double[][] matrix) {
        if (matrix.length == 2) {//if the matrix was 2x2
            return (matrix[0][0] * matrix[1][1]) - (matrix[0][1] * matrix[1][0]);
        } else if (matrix.length == 3) {//if the matrix was 3x3
            return matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
                    - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
                    + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
        }
        return 0;
    }
    //this method is a part of the logic that we do have
    private static double[][] copyMatrix(double[][] matrix) {
        double[][] copy = new double[matrix.length][matrix[0].length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                copy[i][j] = matrix[i][j];
            }
        }
        return copy;
    }


    }

