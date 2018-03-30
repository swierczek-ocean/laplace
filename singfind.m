function [hshift,vshift] = singfind(fun)

%% find the singularities
syms s
P = poles(fun,s);
fun1 = @(s)1/fun(s);
Z = vpasolve(fun1(s) == 0, s);
S = [P;Z];
Q = real(S);
[Q,I] = sort(Q,'descend');
if(isempty(Q)==1)
    hshift = 0;
    vshift = 0;
else
    hshift = abs(Q(1));
    if(size(Q,1)==1)
        vshift = double(imag(S(I(1))));
    elseif(Q(1)>Q(2))
        vshift = double(imag(S(I(1))));
    else
        vshift = double(imag(0.5*(S(I(1))+S(I(2)))));
    end
end
%%

end

