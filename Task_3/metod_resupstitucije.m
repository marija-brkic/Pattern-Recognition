function [s_opt, v0_opt, Neps_opt, M1, M2, S1, S2] = metod_resupstitucije(X1, X2)
    M1 = mean(X1)';
    M2 = mean(X2)';
    S1 = cov(X1);
    S2 = cov(X2);

    s = 0:1e-3:1;
    v0_opt_s = [];
    Neps_s = [];
    
    for k = 1:length(s)
        v = (s(k)*S1 + (1-s(k))*S2)^( -1)*(M2-M1);
        v0 = [];
        Neps = [];
        Y1 = v'*X1';
        Y2 = v'*X2';
        Y = [Y1 Y2]; Y = sort(Y);
        for j = 1:length(Y)-1
            v0(j) = -(Y(j)+Y(j+1))/2;
            Neps(j) = 0;

            for i = 1:length(Y1)
                if (Y1(i) > -v0(j))
                    Neps(j) = Neps(j) +1;
                end
            end
            for i = 1:length(Y2)
                if (Y2(i) < -v0(j))
                    Neps(j) = Neps(j) +1;
                end
            end

        end
        [Neps_s(k), ind] = min(Neps);
        v0_opt_s(k) = v0(ind);

    end
    [Neps_opt, ind] = min(Neps_s);
    v0_opt = v0_opt_s(ind);
    s_opt = s(ind);
end