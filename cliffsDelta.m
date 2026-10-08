function d = cliffsDelta(x, y)
    n = 0;
    for i = 1:length(x)
        for j = 1:length(y)
            n = n + sign(x(i) - y(j));
        end
    end
    d = n / (length(x) * length(y));
end