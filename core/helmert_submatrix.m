function H = helmert_submatrix(K)
H = zeros(K, K-1);
for i=1:(K-1)
    H(1:i, i) =  1 / sqrt(i*(i+1));
    H(i+1, i) = -i / sqrt(i*(i+1));
end
end

