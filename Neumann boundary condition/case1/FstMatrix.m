function SM = FstMatrix(M)

for i = 1 : M
    for j = 1 : M
        SM(i, j) = sin (i * j * pi / (M + 1)) ;
    end
end

SM = sqrt(2) * SM / sqrt(M + 1) ;