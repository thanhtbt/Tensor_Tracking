function er = sub_est_per(U_tr,U_es,type)

    if nargin < 3
        type = 'SEP';
    end
    [n,r] = size(U_tr);
    switch type
        case 'SEP'
            U_es =  orth(U_es);
            U_tr =  orth(U_tr);
            er   =  abs(trace(U_es' * (eye(n)-U_tr*U_tr') * U_es)/trace(U_es'*(U_tr*U_tr')*U_es));
        case 'RE'
            U_tr = orth(U_tr);
            U_es = orth(U_es);
            ER = U_tr - U_es * (U_es'*U_tr); 
            er = norm(ER,'fro');
        case 'Angle'
            er = subspace(U_tr,U_es);  
            er = sin(er);   
        case 'SE' 
            [U_tr, ~] = qr(U_tr);
            [U_es, ~] = qr(U_es); 
            er = norm((eye(n) - U_tr(:, 1 : r)  * U_tr(:, 1 : r)') * U_es(:, 1 : r), 2);
        case 'EV'
            er =  abs(trace(U_tr' * U_es * U_es' * U_tr) /  trace(U_es * U_es')); 
    end

end


