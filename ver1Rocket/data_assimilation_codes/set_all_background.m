function [Xn,XX,Xo,X] = set_all_background(X,Vxi,RXo,timecounterXs) 
%now get the background from model 
Xo = RXo; 
%lets assume that Xn is static for the moment
Xn = Xo;
XX= RXo;
for I = 1:length(Vxi)
Xcnt = timecounterXs+I-1;

RE1 = ['X.Xk',num2str(Xcnt),'=RXo;'];
eval(RE1)

end 
