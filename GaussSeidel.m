%Function to solve a linear system of equations, defined by a square matrix
%A and a vector B, using the Gauss-Seidel method with an initial guess
%defined by an array x_0
%
%The function syntax is x=GaussSeidel(A,B,x_0)
%
%Outputs
%x - n x 1 array - solutions of the linear system, where n is the number of
%                  solutions
%
%Inputs
%A   - n x n array - matrix containing coefficients of the variables for 
%      each equation of the linear system
%B   - n x 1 array - vector containing the constant of each equation of the 
%      system
%x_0 - 1 x n array - array containing initial guesses for each of the 
%      unkowns, used for the Gauss-Seidel method (usually filled with
%      zeros)


function x = GaussSeidel(A,B,x_0)
    
    %check if the provided inputs are of a desired format
    if ~((isnumeric(A) && ismatrix(B)) && (isnumeric(B) && ismatrix(B)) && (isnumeric(x_0) && ismatrix(x_0)))
        error('All inputs of GaussSeidel() must be numeric arrays')
    end
    
    %check if the vector B is of desired shape and size
    if (size(B,2) == size(A,1)) && (size(B,1) == 1)
        B = B.';
    elseif (size(B,1) ~= size(A,1)) || (size(B,2)>1)
        error('The constants array (B) is not correctly dimensioned for a given matrix A')
    end
    
    %check if there are more columns than rows and raise an appropriate 
    %error
    if size(A,1) < size(A,2)
        error('The system of equations has infinitely many solutions')        
    end
    
    
    %create a temporary array A_temp equal to A with B appended on the n+1
    %column
    A_temp = A;
    A_temp(:,size(A,2)+1) = B;
    
    %define the rank of A and A_temp
    rank_A = rank(A);
    rank_A_temp = rank(A_temp);
    
    %define the number of columns of A and assign it to a variable 'n'
    n = size(A,1);
    
    %implement the linear algebra theory connecting rank of matrix A and
    %A_temp to the number of solutions of the linear system
    if rank_A_temp == rank_A && rank_A < n
        if rank_A < size(A,2)
            %display an error according to the theory
            error('The system of equations has infinitely many solutions')
        else
            %if there are more rows than columns, check for duplicate rows:
            
            %define variable temp as the first column of matrix A with 0 elements
            %changed to 1. This will later be used to divide all rows of A and B by
            %the first coefficient of a given row in the linear system
            temp = A(:,1);
            temp(abs(temp) < 10^(-15)) = 1;
            %create a grid of the vector temp repeated over the number of columns
            %of A
            temp = meshgrid(temp,1:size(A,2)).';

            %find the unique rows (ia) of the matrix A divided by the first
            %coefficient in each row, using the unique() function
            [~,ia,~] = unique(A./temp,'rows','stable');
            
            if length(ia) == n
                error('The system of equations has infinitely many solutions')
            else
                %if there exist duplicate rows in the matrix A, inform the
                %user about the possibility of a unique solution after the
                %given rows are deleted
                duplicates = sprintf('%.0f,' , setdiff(1:size(A,1),ia));
                error_message = ['Matrix non-square, but possible unique solution. Duplicate/contradicting rows: ',duplicates(1:end-1)];
                error(error_message)
            end
        end
        
    elseif rank_A < rank_A_temp
        %display an appropriate error, according to theory
        error('System of equations has no solutions')
    end
    
    
    %define the diagonal line of the matrix using diag() function
    diagonal = diag(A);
        
    %find the diagonal elements of A that are equal to zero and assign
    %their position to a variable 'n'
    n = find(abs(diagonal) < 10^(-16));
    
    %if there are any zeros in the diagonal of A, iterate over the rows 
    %with zero diagonal element and try to manipulate the matrix to make
    %the diagonals non-zero 
    for i = n.'         
        available_rows = find(A(:,i)~=0);  %find all rows which have non-zero elements on the position of the current zero diagonal element
        if ~isempty(available_rows)        %check if the array of available rows is not empty
            %change the ith row of A to a sum of the orignal row and first row with a non-zero element in column i 
            A(i,:) = A(i,:) + A(available_rows(1),:);
            B(i,:) = B(i,:) + B(available_rows(1),:);
        else
            error('System of equations cannot be solved using the Gauss-Seidel method') %if the manipulation is not possible Gauss-Seidel method will diverge
        end
    end
    
    %define desired precision of the final answers
    epsilon = 10^(-8);
    
    %define maximum number of iterations as 20000 times the square root of 
    %the size of matrix A - this is due to an assumption that the number of
    %required iterations tends to increase with the square root of the size
    %of the problem
    max_iterations = 20000*sqrt(size(A,2));
    
    %check if array x_0 is of desired format and modify it accordingly
    if (size(x_0,1) == size(A,2)) && (size(x_0,2)==1)
        x_0 = x_0.'; 
    elseif (size(x_0,2) ~= size(A,2)) || (size(x_0,1)>1) %if the array is of an incorrect size, change it to an array of zeros
        x_0 = zeros(1,size(A,2));
    end
    
    %initialize the array 'x' as the guess array x_0
    x = x_0;
    
    %initialize array 'x_1' used to check precision of the answers
    x_1 = x_0 + ones(1,length(x_0));
    
    %make number of iterations initially equal to 0
    iterations = 0;
    
    %iterate while the answers are of insufficient precision 
    %(the abs(x)*epsilon expression makes sure the code isn't susceptible 
    %to small values of x) and while the number of iterations does not 
    %exceed the artificial maximum 
    while any(abs(x - x_1) > abs(x)*epsilon) && (iterations < max_iterations)
        x_1 = x;
        %iterate over rows in A
        for i = 1:size(A,2)
            %for each row estimate x(i) as the negative sum of the 
            %variables (excluding the x(i) variable) multiplied by their 
            %respective coefficients, subtracted from the constant B(i) and
            %divided by the coefficient of x(i) 
            x(i) = (B(i)-sum(A(i,[1:i-1,i+1:size(A,2)]).*x([1:i-1,i+1:size(x,2)])))/A(i,i);
        end
        
        %check if any of the computed solutions are non-numeric or infinite
        %to stop the program from running unnecessarily
        if any(isnan(x)) || any(isinf(x))
            error('Linear system diverges under the Gauss-Seidel method')
        end
        
        iterations = iterations+1; %increment the iterations variable
    end
            
    %assign 1 (true) to the variable diagonally_dominant
    diagonally_dominant = 1;
    
    %check if matrix A is diagonally dominant (the absolute value of the
    %diagonal entry is larger than or equal to the sum of the other entries
    if any(abs(diag(A)) <= (sum(abs(A),2)-diag(A)))
        diagonally_dominant = 0;
    end
    
    %use the chol() function to find wheter the matrix A is positive
    %definite - p equals to 0 when positive definite and p>0 otherwise
    [~,p] = chol(A);
    
    %check if the method has reached the maximum number of iterations
    if iterations >= max_iterations
        %if satisfied we suspect that the matrix might have diverged, so we
        %check if matrix is not positive definite or not diagonally 
        %dominant (conditions opposite to ones guaranteeing Gauss-Seidel 
        %convergence)
        if p > 0 || diagonally_dominant == 0
            error('Linear system diverges or iterates in a circle') 
        else
            %otherwise display a warning
            warning('Maximum number of iterations reached - solutions might be inaccurate')
        end
    end