% given an array of identity vectors [n,m] where m columns indicate unique
% subjects and n rows are the associated observations. compute the cosine
% distance between each possible pairing of i-Vectors. the resultant matrix
% will be of size m*m which half of it being redudant. the formula for
% computing the cosine distance comes from a talk given by Howard Lei on
% Joint Factor Analysis (JFA) and i-vector Tutorial on behalf of the AFRL
% ICSI. if the ivectors point in the same direction the result should be 1
% and if they point in opposite directions the value should be -1. 1 and -1
% being the highest and lowest values possible.
function result = cosineDistance(varargin)
a0 = tic;
if( nargin == 1 )
    % evaluate EVERYTHING
    i_vector = varargin{1};
    [~,columns] = size(i_vector);
    result = ones(columns,columns);
    
    parfor i=1:columns
        target = i_vector(:,i);
        target_norm = norm(target);
        for k=1:columns
            result(i,k) = ( target'*i_vector(:,k) ) / ...
                ( target_norm*norm(i_vector(:,k)) );
        end
    end
elseif( nargin == 2)
    % only produce one quadrant, testing INTO training
    enroll_iv = varargin{1};
    test_iv = varargin{2};
    [~, n_enroll] = size(enroll_iv);
    [~, n_test] = size(test_iv);
    target_norm = zeros(n_test,1);
    enroll_norm = zeros(n_enroll,1);
    for i=1:n_test
        target_norm(i) = norm(test_iv(:,i));
    end
    for i=1:n_enroll
        enroll_norm(i) = norm(enroll_iv(:,i));
    end
    result = zeros(n_enroll,n_test);
    for t=1:n_test
        target = test_iv(:,t);
        for e=1:n_enroll
            result(e,t) = ( target'*enroll_iv(:,e) ) / ...
                ( target_norm(t) * enroll_norm(e) );
        end
    end
end
a1 = toc(a0);
% fprintf('Cosine Distance completed in %f seconds.\n', a1);
end