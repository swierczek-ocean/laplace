function F = master_laplace_fcn(s,a,b,jj)
%% description
% index jj selects the jjth laplace function.
% the function is evaluated at complex frequencies s
% with parameter a and possible b
%%

if(jj==1)
    F = 1./(s-a);
elseif(jj==2)
    F = 1./(s-a).^2;
elseif(jj==3)
    F = 1./(s-a).^3;
elseif(jj==4)
    F = 1./(s-a).^4;
elseif(jj==5)
    F = 1./(s-a).^5;
elseif(jj==6)
    F = 1./(s.^2 +a^2);
elseif(jj==7)
    F = s./(s.^2 +a^2);
elseif(jj==8)
    F = 1./(s.^2 - a^2);
elseif(jj==9)
    F = s./(s.^2 - a^2);
elseif(jj==10)
    F = 1./(s.^2 +a^2).^2;
elseif(jj==11)
    F = s./(s.^2 +a^2).^2;
elseif(jj==12)
    F = (s.^2)./(s.^2 +a^2).^2;
elseif(jj==13)
    F = (s.^3)./(s.^2 +a^2).^2;
elseif(jj==14)
    F = (s.^2-a^2)./(s.^2 +a^2).^2;
elseif(jj==15)
    F = 1./(s.^2 - a^2).^2;
elseif(jj==16)
    F = s./(s.^2 - a^2).^2;
elseif(jj==17)
    F = s.^2./(s.^2 - a^2).^2;
elseif(jj==18)
    F = s.^3./(s.^2 - a^2).^2;
elseif(jj==19)
    F = (s.^2+a^2)./(s.^2 - a^2).^2;
elseif(jj==20)
    F = 1./(s.^2+a^2).^3;
elseif(jj==21)
    F = s./(s.^2+a^2).^3;
elseif(jj==22)
    F = s.^2./(s.^2+a^2).^3;
elseif(jj==23)
    F = s.^3./(s.^2+a^2).^3;
elseif(jj==24)
    F = s.^4./(s.^2+a^2).^3;
elseif(jj==25)
    F = s.^5./(s.^2+a^2).^3;
elseif(jj==26)
    F = (3.*s.^2-a^2)./(s.^2+a^2).^3;
elseif(jj==27)
    F = (s.^3-3*a^2.*s)./(s.^2+a^2).^3;
elseif(jj==28)
    F = (s.^4-6*a^2.*s.^2+a^4)./(s.^2+a^2).^4;
elseif(jj==29)
    F = (s.^3 - a^2.*s)./(s.^2+a^2).^4;
elseif(jj==30)
    F = 1./(s.^2-a^2).^3;
elseif(jj==31)
    F = s./(s.^2-a^2).^3;
elseif(jj==32)
    F = s.^2./(s.^2-a^2).^3;
elseif(jj==33)
    F = s.^3./(s.^2-a^2).^3;
elseif(jj==34)
    F = s.^4./(s.^2-a^2).^3;
elseif(jj==35)
    F = s.^5./(s.^2-a^2).^3;
elseif(jj==36)
    F = (3.*s.^2+a^2)./(s.^2-a^2).^3;
elseif(jj==37)
    F = (s.^3+3*a^2.*s)./(s.^2-a^2).^3;
elseif(jj==38)
    F = (s.^4 +6*a^2.*s.^2+a^4)./(s.^2-a^2).^4;
elseif(jj==39)
    F = (s.^3+a^2.*s)./(s.^2-a^2).^4;
elseif(jj==40)
    F = 1./(s.^3+a^3);
elseif(jj==41)
    F = s./(s.^3+a^3);
elseif(jj==42)
    F = s.^2./(s.^3+a^3);
elseif(jj==43)
    F = 1./(s.^3-a^3);
elseif(jj==44)
    F = s./(s.^3-a^3);
elseif(jj==45)
    F = s.^2./(s.^3-a^3);
elseif(jj==46)
    F = 1./(s.^4+a^4);
elseif(jj==47)
    F = s./(s.^4+a^4);
elseif(jj==48)
    F = s.^2./(s.^4+a^4);
elseif(jj==49)
    F = s.^3./(s.^4+a^4);
elseif(jj==50)
    F = 1./(s.^4-a^4);
elseif(jj==51)
    F = s./(s.^4-a^4);
elseif(jj==52)
    F = s.^2./(s.^4-a^4);
elseif(jj==53)
    F = s.^3./(s.^4-a^4);
elseif(jj==54)
    F = 1./(s.*sqrt(s+a));
elseif(jj==55)
    F = 1./(sqrt(s).*(s-a));
elseif(jj==56)
    F = 1./sqrt(s.^2+a^2);
elseif(jj==57)
    F = 1./sqrt(s.^2-a^2);
elseif(jj==58)
    F = (sqrt(s.^2+a^2)-s)./sqrt(s.^2+a^2);
elseif(jj==59)
    F = (sqrt(s.^2+a^2)-s).^2./sqrt(s.^2+a^2);
elseif(jj==60)
    F = (sqrt(s.^2+a^2)-s).^3./sqrt(s.^2+a^2);
elseif(jj==61)
    F = (sqrt(s.^2+a^2)-s).^4./sqrt(s.^2+a^2);
elseif(jj==62)
    F = (s - sqrt(s.^2-a^2))./sqrt(s.^2-a^2);
elseif(jj==63)
    F = (s - sqrt(s.^2-a^2)).^2./sqrt(s.^2-a^2);
elseif(jj==64)
    F = (s - sqrt(s.^2-a^2)).^3./sqrt(s.^2-a^2);
elseif(jj==65)
    F = (s - sqrt(s.^2-a^2)).^4./sqrt(s.^2-a^2);
elseif(jj==66)
    F = 1./(s.^2+a^2).^(3/2);
elseif(jj==67)
    F = s./(s.^2+a^2).^(3/2);
elseif(jj==68)
    F = s.^2./(s.^2+a^2).^(3/2);
elseif(jj==69)
    F = 1./(s.^2-a^2).^(3/2);
elseif(jj==70)
    F = s./(s.^2-a^2).^(3/2);
elseif(jj==71)
    F = s.^2./(s.^2-a^2).^(3/2);
elseif(jj==72)
    F = 1./(s.*(exp(s)-a));
elseif(jj==73)
    F = (exp(s)-1)./(s.*(exp(s)-a));
elseif(jj==74)
    F = exp(-a./s)./sqrt(s);
elseif(jj==75)
    F = exp(-a./s)./(s.^(3/2));
elseif(jj==76)
    F = exp(-a./s)./s;
elseif(jj==77)
    F = exp(-a./s)./(s.^2);
elseif(jj==78)
    F = exp(-a./s)./(s.^3);
elseif(jj==79)
    F = exp(-a./s)./(s.^4);
elseif(jj==80)
    F = exp(-a.*sqrt(s))./sqrt(s);
elseif(jj==81)
    F = exp(-a.*sqrt(s));
elseif(jj==82)
    F = (1-exp(-a.*sqrt(s)))./s;
elseif(jj==83)
    F = exp(-a.*sqrt(s))./s;
elseif(jj==84)
    F = exp(-a./sqrt(s))./s;
elseif(jj==85)
    F = exp(-a./sqrt(s))./(s.^2);
elseif(jj==86)
    F = exp(-a./sqrt(s))./(s.^3);
elseif(jj==87)
    F = exp(-a./sqrt(s))./(s.^4);
elseif(jj==88)
    F = exp(-a./sqrt(s))./(s.^5);
elseif(jj==89)
    F = 0.5.*log((s.^2+a^2)/(a^2))./s;
elseif(jj==90)
    F = log((s+a)./a)./s;
elseif(jj==91)
    
elseif(jj==92)
    
elseif(jj==93)
    
elseif(jj==94)

elseif(jj==95)
    
elseif(jj==96)
    
elseif(jj==97)
    
elseif(jj==98)
    
elseif(jj==99)
    
elseif(jj==100)
    
elseif(jj==101)
    
elseif(jj==102)
    
elseif(jj==103)
    
elseif(jj==104)
    
elseif(jj==105)

elseif(jj==106)
    
elseif(jj==107)
    
elseif(jj==108)
    
elseif(jj==109)
    
elseif(jj==110)
    
elseif(jj==111)
    
elseif(jj==112)
    
elseif(jj==113)
    
elseif(jj==114)

elseif(jj==115)
    
elseif(jj==116)
    
elseif(jj==117)
    
elseif(jj==118)
    
elseif(jj==119)
    
elseif(jj==120)
    
elseif(jj==121)
    
elseif(jj==122)
    
elseif(jj==123)
    
elseif(jj==124)
    
elseif(jj==125)

elseif(jj==126)
    
elseif(jj==127)
    
elseif(jj==128)
    
elseif(jj==129)
    
elseif(jj==130)
    
elseif(jj==131)
    
elseif(jj==132)
    
elseif(jj==133)
    
elseif(jj==134)

elseif(jj==135)
    
elseif(jj==136)
    
elseif(jj==137)
    
elseif(jj==138)
    
elseif(jj==139)
    
elseif(jj==140)
    
elseif(jj==141)
    
elseif(jj==142)
    
elseif(jj==143)
    
elseif(jj==144)
    
elseif(jj==145)

elseif(jj==146)
    
elseif(jj==147)
    
elseif(jj==148)
    
elseif(jj==149)
    
elseif(jj==150)
    
elseif(jj==151)
    
elseif(jj==152)
    
elseif(jj==153)
    
elseif(jj==154)

elseif(jj==155)
    
elseif(jj==156)
    
elseif(jj==157)
    
elseif(jj==158)
    
elseif(jj==159)
    
elseif(jj==160)
    
elseif(jj==161)
    
elseif(jj==162)
    
elseif(jj==163)
    
elseif(jj==164)
    
elseif(jj==165)

elseif(jj==166)
    
elseif(jj==167)
    
elseif(jj==168)
    
elseif(jj==169)
    
elseif(jj==170)
    
elseif(jj==171)
    
elseif(jj==172)
    
elseif(jj==173)
    
elseif(jj==174)

elseif(jj==175)
    
elseif(jj==176)
    
elseif(jj==177)
    
elseif(jj==178)
    
elseif(jj==179)
    
elseif(jj==180)
    
elseif(jj==181)
    
elseif(jj==182)
    
elseif(jj==183)
    
elseif(jj==184)
    
elseif(jj==185)

elseif(jj==186)
    
elseif(jj==187)
    
elseif(jj==188)
    
elseif(jj==189)
    
elseif(jj==190)
    
elseif(jj==191)
    
elseif(jj==192)
    
elseif(jj==193)
    
elseif(jj==194)
    
elseif(jj==195)

elseif(jj==196)
    
elseif(jj==197)
    
elseif(jj==198)
    
elseif(jj==199)
    
elseif(jj==200)  

    

end

