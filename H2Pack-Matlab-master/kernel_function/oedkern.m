function KXX = oedkern(X_input, ls, sig_f)

if isnumeric(X_input)
    KXX = gaussKern(X_input, sig_f, ls);
elseif iscell(X_input)&&length(X_input)==2     
    KXX = gaussKern(X_input{1}, sig_f, ls, X_input{2});
end
end