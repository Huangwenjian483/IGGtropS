function S = extractCoefficients(pathB, pathV, ihgt, ilat)
%==========================input============================
% pathB: The file path for storing the mask matrix B
% pathV: The file path for storing the non-zero coefficients V
% ihgt: height layer index of target point
% ilat: latitude band index of target point
%===========================================================
%==========================output===========================
% S: The required sparse model coefficients (144*4)
%===========================================================
    n_row = 144;
    page_size = 144 * 73;

    idx1 = (ihgt - 1) * page_size + (ilat - 1) * n_row + (1 : 2 * n_row);
    idx2 = ihgt       * page_size + (ilat - 1) * n_row + (1 : 2 * n_row);
    maxidx = idx2(end);

    fid = fopen(pathB, 'r');
    B = fread(fid, maxidx, 'ubit1');
    fclose(fid);

    B = logical(B);
    B1 = B(idx1);
    B2 = B(idx2);

    sum1 = 0;
    if ~(ihgt == 1 && ilat == 1)
        sum1 = sum(B(1 : idx1(1) - 1));
    end
    count1 = sum(B1);
    sum2 = sum1 + count1 + sum(B(idx1(end) + 1 : idx2(1) - 1));
    count2 = sum(B2);

    fid = fopen(pathV);
    V = fread(fid, sum2+count2, 'float');
    fclose(fid);

    V1 = V(sum1 + 1 : sum1 + count1);
    V2 = V(sum2 + 1 : sum2 + count2);

    S1 = zeros(2 * n_row, 1); 
    S1(B1) = V1;
    S2 = zeros(2 * n_row, 1); 
    S2(B2) = V2;

    S = reshape([S1; S2], n_row, 4);
end