%=========================================================================
%                   Discretize the Income Process
%=========================================================================
% Description: Builds the Lambda matrix associated with the jump-drift
%              income process
% 
%=========================================================================


function Lambda_matrix = Build_Lambda_matrix(lambda_z,n_z,n_a)

sparse_lambda_z = [];
for m=1:n_z % row
    row = [];
    for j=1:n_z % column
        block = lambda_z(m,j)*speye(n_a);
        row = [row,block];
    end    
    sparse_lambda_z = vertcat(sparse_lambda_z,row);
end
Lambda_matrix = sparse_lambda_z;

end
