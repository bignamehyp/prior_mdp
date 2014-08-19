function sse=fitDDM(params,x,y)
A=params(1);
k=params(2);
Fitted_Curve= A/k./x .* tanh(A*k*x);
Error_Vector=Fitted_Curve - y;
% When curvefitting, a typical quantity to
% minimize is the sum of squares error
sse=sum(Error_Vector.^2);
% You could also write sse as
% sse=Error_Vector(:)'*Error_Vector(:);