function Sings = makesings(a,b)
Sings = zeros(200,2);
Sings(1:5,1) = a*ones(5,1);
Sings(8:9,1) = a*ones(2,1);
Sings(15:19,1) = a*ones(5,1);
Sings(30:39,1) = a*ones(10,1);
Sings(40:42,1) = 0.5*a*ones(3,1);
Sings(43:53,1) = a*ones(11,1);
Sings(55,1) = a;
Sings(57,1) = a;
Sings(62:65,1) = a*ones(4,1);
Sings(69:71,1) = a*ones(3,1);
Sings(72:73,1) = log(a)*ones(2,1);
Sings(122,1) = log(a);
Sings(125,1) = a;
Sings(128,1) = -b;
Sings(130:131,1) = b*ones(2,1);
Sings(132:133,1) = (b+a)*ones(2,1);
Sings(134:135,1) = max([a,b])*ones(2,1);












end

