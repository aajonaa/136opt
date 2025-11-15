
function [lb,ub,dim,fobj] = Get_Functions(F,DimValue)
if nargin == 1
    DimValue = 10;
end

switch F
    case 'F1'
        fobj = @F1;
        lb=-100;
        ub=100;
        dim=DimValue;
        
    case 'F2'
        fobj = @F2;
        lb=-10;
        ub=10;
        dim=DimValue;
        
    case 'F3'
        fobj = @F3;
        lb=-100;
        ub=100;
        dim=DimValue;
        
    case 'F4'
        fobj = @F4;
        lb=-100;
        ub=100;
        dim=DimValue;
        
    case 'F5'
        fobj = @F5;
        lb=-30;
        ub=30;
        dim=DimValue;
        
    case 'F6'
        fobj = @F6;
        lb=-100;
        ub=100;
        dim=DimValue;
        
    case 'F7'
        fobj = @F7;
        lb=-1.28;
        ub=1.28;
        dim=DimValue;
        
    case 'F8'
        fobj = @F8;
        lb=-500;
        ub=500;
        dim=DimValue;
        
    case 'F9'
        fobj = @F9;
        lb=-5.12;
        ub=5.12;
        dim=DimValue;
        
    case 'F10'
        fobj = @F10;
        lb=-32;
        ub=32;
        dim=DimValue;
        
    case 'F11'
        fobj = @F11;
        lb=-600;
        ub=600;
        dim=DimValue;
        
    case 'F12'
        fobj = @F12;
        lb=-50;
        ub=50;
        dim=DimValue;
        
    case 'F13'
        fobj = @F13;
        lb=-50;
        ub=50;
        dim=DimValue;
        
    case 'F14'
        fobj = @F14;
        lb=-65.536;
        ub=65.536;
        dim=2;
        
    case 'F15'
        fobj = @F15;
        lb=-5;
        ub=5;
        dim=4;
        
    case 'F16'
        fobj = @F16;
        lb=-5;
        ub=5;
        dim=2;
        
    case 'F17'
        fobj = @F17;
        lb=-5;
        ub=5;
        dim=2;
        
    case 'F18'
        fobj = @F18;
        lb=-2;
        ub=2;
        dim=2;
        
    case 'F19'
        fobj = @F19;
        lb=0;
        ub=1;
        dim=3;
        
    case 'F20'
        fobj = @F20;
        lb=0;
        ub=1;
        dim=6;     
        
    case 'F21'
        fobj = @F21;
        lb=0;
        ub=10;
        dim=4;    
        
    case 'F22'
        fobj = @F22;
        lb=0;
        ub=10;
        dim=4;    
        
    case 'F23'
        fobj = @F23;
        lb=0;
        ub=10;
        dim=4;  
    case 'F24'
        fobj = @F24;
        lb=-100;
        ub=100;
        dim=DimValue;  %Test functions are only defined for D=2,10,30,50,100.
    case 'F25'
        fobj = @F25;
        lb=-100;
        ub=100;
        dim=DimValue;  %Test functions are only defined for D=2,10,30,50,100.
    case 'F26'
        fobj = @F26;
        lb=-100;
        ub=100;
        dim=DimValue;  %Test functions are only defined for D=2,10,30,50,100.
    case 'F27'
        fobj = @F27;
        lb=-100;
        ub=100;
        dim=DimValue;  %Test functions are only defined for D=2,10,30,50,100.
    case 'F28'
        fobj = @F28;
        lb=-100;
        ub=100;
        dim=DimValue;  %Test functions are only defined for D=2,10,30,50,100.
    case 'F29'
        fobj = @F29;
        lb=-100;
        ub=100;
        dim=DimValue;  %Test functions are only defined for D=2,10,30,50,100.
    case 'F30'
        fobj = @F30;
        lb=-100;
        ub=100;
        dim=DimValue;
    case 'F31'
        fobj = @F31;
        lb=-100;
        ub=100;
        dim=DimValue;  %Test functions are only defined for D=2,10,30,50,100.
    case 'F32'
        fobj = @F32;
        lb=-100;
        ub=100;
        dim=DimValue;  %Test functions are only defined for D=2,10,30,50,100.
    case 'F33'
        fobj = @F33;
        lb=-100;
        ub=100;
        dim=DimValue;  %Test functions are only defined for D=2,10,30,50,100.
    case 'F34'
        fobj = @F34;
        lb=-100;
        ub=100;
        dim=DimValue;  %Test functions are only defined for D=2,10,30,50,100.
    case 'F35'
        fobj = @F35;
        lb=-100;
        ub=100;
        dim=DimValue;  %Test functions are only defined for D=2,10,30,50,100.
    case 'F36'
        fobj = @F36;
        lb=-100;
        ub=100;
        dim=DimValue;  %Test functions are only defined for D=2,10,30,50,100.
    case 'F37'
        fobj = @F37;
        lb=-100;
        ub=100;
        dim=DimValue;  %Test functions are only defined for D=2,10,30,50,100.
    case 'F38'
        fobj = @F38;
        lb=-100;
        ub=100;
        dim=DimValue;  %Test functions are only defined for D=2,10,30,50,100.
    case 'F39'
        fobj = @F39;
        lb=-100;
        ub=100;
        dim=DimValue;  %Test functions are only defined for D=2,10,30,50,100.
    case 'F40'
        fobj = @F40;
        lb=-100;
        ub=100;
        dim=DimValue;
    case 'F41'
        fobj = @F41;
        lb=-100;
        ub=100;
        dim=DimValue;  %Test functions are only defined for D=2,10,30,50,100.
    case 'F42'
        fobj = @F42;
        lb=-100;
        ub=100;
        dim=DimValue;  %Test functions are only defined for D=2,10,30,50,100.
    case 'F43'
        fobj = @F43;
        lb=-100;
        ub=100;
        dim=DimValue;  %Test functions are only defined for D=2,10,30,50,100.
    case 'F44'
        fobj = @F44;
        lb=-100;
        ub=100;
        dim=DimValue;  %Test functions are only defined for D=2,10,30,50,100.
    case 'F45'
        fobj = @F45;
        lb=-100;
        ub=100;
        dim=DimValue;  %Test functions are only defined for D=2,10,30,50,100.
    case 'F46'
        fobj = @F46;
        lb=-100;
        ub=100;
        dim=DimValue;  %Test functions are only defined for D=2,10,30,50,100.
    case 'F47'
        fobj = @F47;
        lb=-100;
        ub=100;
        dim=DimValue;  %Test functions are only defined for D=2,10,30,50,100.
    case 'F48'
        fobj = @F48;
        lb=-100;
        ub=100;
        dim=DimValue;  %Test functions are only defined for D=2,10,30,50,100.
    case 'F49'
        fobj = @F49;
        lb=-100;
        ub=100;
        dim=DimValue;  %Test functions are only defined for D=2,10,30,50,100.
    case 'F50'
        fobj = @F50;
        lb=-100;
        ub=100;
        dim=DimValue;
    case 'F51'
        fobj = @F51;
        lb=-100;
        ub=100;
        dim=DimValue;  %Test functions are only defined for D=2,10,30,50,100.
    case 'F52'
        fobj = @F52;
        lb=-100;
        ub=100;
        dim=DimValue;  %Test functions are only defined for D=2,10,30,50,100.
    case 'F53'
        fobj = @F53;
        lb=-100;
        ub=100;
        dim=DimValue;  %Test functions are only defined for D=2,10,30,50,100.   
    case 'F54'
        fobj = @F54;
        lb=-100;
        ub=100;
        dim=DimValue; 
    case 'F55'
        fobj = @F55;
        lb=-100;
        ub=100;
        dim=DimValue; 
    case 'F56'
        fobj = @F56;
        lb=-100;
        ub=100;
        dim=DimValue;
    case 'F57'
        fobj = @F57;
        lb=-100;
        ub=100;
        dim=DimValue; 
    case 'F58'
        fobj = @F58;
        lb=-100;
        ub=100;
        dim=DimValue; 
    case 'F59'
        fobj = @F59;
        lb=-100;
        ub=100;
        dim=DimValue; 
    case 'F60'
        fobj = @F60;
        lb=0;
        ub=600;
        dim=DimValue;
    case 'F61'
        fobj = @F61;
        lb=-32;
        ub=32;
        dim=DimValue; 
    case 'F62'
        fobj = @F62;
        lb=-5;
        ub=5;
        dim=DimValue;  
    case 'F63'
        fobj = @F63;
        lb=-5;
        ub=5;
        dim=DimValue; 
    case 'F64'
        fobj = @F64;
        lb=-0.5;
        ub=0.5;
        dim=DimValue;  
    case 'F65'
        fobj = @F65;
        lb=-pi;
        ub=pi;
        dim=DimValue; 
    case 'F66'
        fobj = @F66;
        lb=-3;
        ub=1;
        dim=DimValue;  
    case 'F67'
        fobj = @F67;
        lb=-100;
        ub=100;
        dim=DimValue; 
    case 'F68'
        fobj = @F68;
        lb=-5;
        ub=5;
        dim=DimValue;  
    case 'F69'
        fobj = @F69;
        lb=-5;
        ub=5;
        dim=DimValue; 
    case 'F70'
        fobj = @F70;
        lb=-5;
        ub=5;
        dim=DimValue;
    case 'F71'
        fobj = @F71;
        lb=-5;
        ub=5;
        dim=DimValue;  
    case 'F72'
        fobj = @F72;
        lb=-5;
        ub=5;
        dim=DimValue;  
    case 'F73'
        fobj = @F73;
        lb=-5;
        ub=5;
        dim=DimValue;  
    case 'F74'
        fobj = @F74;
        lb=-5;
        ub=5;
        dim=DimValue;  
    case 'F75'
        fobj = @F75;
        lb=-5;
        ub=5;
        dim=DimValue;  
    case 'F76'
        fobj = @F76;
        lb=-5;
        ub=5;
        dim=DimValue;  
    case 'F77'
        fobj = @F77;
        lb=-5;
        ub=5;
        dim=DimValue;  
    case 'F78'
        fobj = @F78;
        lb=-2;
        ub=5;
        dim=DimValue;  
    case 'F79'
        fobj = @F79;
        lb=-100;
        ub=100;
        dim=DimValue;  
    case 'F80'
        fobj = @F80;
        lb=-100;
        ub=100;
        dim=DimValue;  
    case 'F81'
        fobj = @F81;
        lb=-100;
        ub=100;
        dim=DimValue;  
    case 'F82'
        fobj = @F82;
        lb=-100;
        ub=100;
        dim=DimValue; 
    case 'F83'
        fobj = @F83;
        lb=-100;
        ub=100;
        dim=DimValue;  
    case 'F84'
        fobj = @F84;
        lb=-100;
        ub=100;
        dim=DimValue; 
    case 'F85'
        fobj = @F85;
        lb=-100;
        ub=100;
        dim=DimValue;
    case 'F86'
        fobj = @F86;
        lb=-100;
        ub=100;
        dim=DimValue;
    case 'F87'
        fobj = @F87;
        lb=-100;
        ub=100;
        dim=DimValue;  
    case 'F88'
        fobj = @F88;
        lb=-100;
        ub=100;
        dim=DimValue;  
    case 'F89'
        fobj = @F89;
        lb=-100;
        ub=100;
        dim=DimValue; 
    case 'F90'
        fobj = @F90;
        lb=-100;
        ub=100;
        dim=DimValue; 
    case 'F91'
        fobj = @F91;
        lb=-100;
        ub=100;
        dim=DimValue; 
    case 'F92'
        fobj = @F92;
        lb=-100;
        ub=100;
        dim=DimValue;  
    case 'F93'
        fobj = @F93;
        lb=-100;
        ub=100;
        dim=DimValue;  
    case 'F94'
        fobj = @F94;
        lb=-100;
        ub=100;
        dim=DimValue;  
    case 'F95'
        fobj = @F95;
        lb=-100;
        ub=100;
        dim=DimValue;
    case 'F96'
        fobj = @F96;
        lb=-100;
        ub=100;
        dim=DimValue;  
    case 'F97'
        fobj = @F97;
        lb=-100;
        ub=100;
        dim=DimValue; 
    case 'F98'
        fobj = @F98;
        lb=-100;
        ub=100;
        dim=DimValue; 
    case 'F99'
        fobj = @F99;
        lb=-100;
        ub=100;
        dim=DimValue;  
    case 'F100'
        fobj = @F100;
        lb=-100;
        ub=100;
        dim=DimValue;  
    case 'F101'
        fobj = @F101;
        lb=-100;
        ub=100;
        dim=DimValue;  
    case 'F102'
        fobj = @F102;
        lb=-100;
        ub=100;
        dim=DimValue;  
    case 'F103'
        fobj = @F103;
        lb=-100;
        ub=100;
        dim=DimValue;  
    case 'F104'
        fobj = @F104;
        lb=-100;
        ub=100;
        dim=DimValue;  
    case 'F105'
        fobj = @F105;
        lb=-100;
        ub=100;
        dim=DimValue;
    case 'F106'
        fobj = @F106;
        lb=-100;
        ub=100;
        dim=DimValue;  
    case 'F107'
        fobj = @F107;
        lb=-100;
        ub=100;
        dim=DimValue; 
    case 'F108'
        fobj = @F108;
        lb=-100;
        ub=100;
        dim=DimValue; 
    case 'F109'
        fobj = @F109;
        lb=-100;
        ub=100;
        dim=DimValue; 
    case 'F110'
        fobj = @F110;
        lb=-100;
        ub=100;
        dim=DimValue; 
    case 'F111'
        fobj = @F111;
        lb=-100;
        ub=100;
        dim=DimValue;
    case 'F112'
        fobj = @F112;
        lb=-100;
        ub=100;
        dim=DimValue; 
    case 'F113'
        fobj = @F113;
        lb=-100;
        ub=100;
        dim=DimValue; 
    case 'F114'
        fobj = @F114;
        lb=-100;
        ub=100;
        dim=DimValue; 
    case 'F115'
        fobj = @F115;
        lb=-100;
        ub=100;
        dim=DimValue; 
    case 'F116'
        fobj = @F116;
        lb=-100;
        ub=100;
        dim=DimValue; 
    case 'F117'
        fobj = @F117;
        lb=-100;
        ub=100;
        dim=DimValue; 
    case 'F118'
        fobj = @F118;
        lb=-100;
        ub=100;
        dim=DimValue; 
    case 'F119'
        fobj = @F119;
        lb=-100;
        ub=100;
        dim=DimValue; 
    case 'F120'
        fobj = @F120;
        lb=-100;
        ub=100;
        dim=DimValue; 
    case 'F121'
        fobj = @F121;
        lb=-100;
        ub=100;
        dim=DimValue; 
    case 'F122'
        fobj = @F122;
        lb=-100;
        ub=100;
        dim=DimValue; 
    case 'F123'
        fobj = @F123;
        lb=-100;
        ub=100;
        dim=DimValue; 
    case 'F124'
        fobj = @F124;
        lb=-100;
        ub=100;
        dim=DimValue; 
    case 'F125'
        fobj = @F125;
        lb=-100;
        ub=100;
        dim=DimValue; 
    case 'F126'
        fobj = @F126;
        lb=-100;
        ub=100;
        dim=DimValue; 
    case 'F127'
        fobj = @F127;
        lb=-100;
        ub=100;
        dim=DimValue; 
    case 'F128'
        fobj = @F128;
        lb=-100;
        ub=100;
        dim=DimValue; 
    case 'F129'
        fobj = @F129;
        lb=-100;
        ub=100;
        dim=DimValue; 
    case 'F130'
        fobj = @F130;
        lb=-100;
        ub=100;
        dim=DimValue; 
    case 'F131'
        fobj = @F131;
        lb=-100;
        ub=100;
        dim=DimValue; 
    case 'F132'
        fobj = @F132;
        lb=-100;
        ub=100;
        dim=DimValue; 
    case 'F133'
        fobj = @F133;
        lb=-100;
        ub=100;
        dim=DimValue; 
    case 'F134'
        fobj = @F134;
        lb=-100;
        ub=100;
        dim=DimValue; 
    case 'F135'
        fobj = @F135;
        lb=-100;
        ub=100;
        dim=DimValue; 
    case 'F136'
        fobj = @F136;
        lb=-100;
        ub=100;
        dim=DimValue; 
    case 'F137'
        fobj = @F137;
        lb=-8192;
        ub=8192;
        dim=9; 
    case 'F138'
        fobj = @F138;
        lb=-16384;
        ub=16384;
        dim=16; 
    case 'F139'
        fobj = @F139;
        lb=-4;
        ub=4;
        dim=18; 
    case 'F140'
        fobj = @F140;
        lb=-100;
        ub=100;
        dim=10; 
    case 'F141'
        fobj = @F141;
        lb=-100;
        ub=100;
        dim=10; 
    case 'F142'
        fobj = @F142;
        lb=-100;
        ub=100;
        dim=10; 
    case 'F143'
        fobj = @F143;
        lb=-100;
        ub=100;
        dim=10; 
    case 'F144'
        fobj = @F144;
        lb=-100;
        ub=100;
        dim=10; 
    case 'F145'
        fobj = @F145;
        lb=-100;
        ub=100;
        dim=10; 
    case 'F146'
        fobj = @F146;
        lb=-100;
        ub=100;
        dim=10; 
end

end

% F1

function o = F1(x)
o=sum(x.^2);
end

% F2

function o = F2(x)
o=sum(abs(x))+prod(abs(x));
end

% F3

function o = F3(x)
dim=size(x,2);
o=0;
for i=1:dim
    o=o+sum(x(1:i))^2;
end
end

% F4

function o = F4(x)
o=max(abs(x));
end

% F5

function o = F5(x)
dim=size(x,2);
o=sum(100*(x(2:dim)-(x(1:dim-1).^2)).^2+(x(1:dim-1)-1).^2);
end

% F6

function o = F6(x)
o=sum(abs((x+.5)).^2);
end

% F7

function o = F7(x)
dim=size(x,2);
o=sum([1:dim].*(x.^4))+rand;
end

% F8

function o = F8(x)
o=sum(-x.*sin(sqrt(abs(x))));
end

% F9

function o = F9(x)
dim=size(x,2);
o=sum(x.^2-10*cos(2*pi.*x))+10*dim;
end

% F10

function o = F10(x)
dim=size(x,2);
o=-20*exp(-.2*sqrt(sum(x.^2)/dim))-exp(sum(cos(2*pi.*x))/dim)+20+exp(1);
end

% F11

function o = F11(x)
dim=size(x,2);
o=sum(x.^2)/4000-prod(cos(x./sqrt([1:dim])))+1;
end

% F12

function o = F12(x)
dim=size(x,2);
o=(pi/dim)*(10*((sin(pi*(1+(x(1)+1)/4)))^2)+sum((((x(1:dim-1)+1)./4).^2).*...
(1+10.*((sin(pi.*(1+(x(2:dim)+1)./4)))).^2))+((x(dim)+1)/4)^2)+sum(Ufun(x,10,100,4));
end

% F13

function o = F13(x)
dim=size(x,2);
o=.1*((sin(3*pi*x(1)))^2+sum((x(1:dim-1)-1).^2.*(1+(sin(3.*pi.*x(2:dim))).^2))+...
((x(dim)-1)^2)*(1+(sin(2*pi*x(dim)))^2))+sum(Ufun(x,5,100,4));
end

% F14

function o = F14(x)
aS=[-32 -16 0 16 32 -32 -16 0 16 32 -32 -16 0 16 32 -32 -16 0 16 32 -32 -16 0 16 32;,...
-32 -32 -32 -32 -32 -16 -16 -16 -16 -16 0 0 0 0 0 16 16 16 16 16 32 32 32 32 32];

for j=1:25
    bS(j)=sum((x'-aS(:,j)).^6);
end
o=(1/500+sum(1./([1:25]+bS))).^(-1);
end

% F15

function o = F15(x)
aK=[.1957 .1947 .1735 .16 .0844 .0627 .0456 .0342 .0323 .0235 .0246];
bK=[.25 .5 1 2 4 6 8 10 12 14 16];bK=1./bK;
o=sum((aK-((x(1).*(bK.^2+x(2).*bK))./(bK.^2+x(3).*bK+x(4)))).^2);
end

% F16

function o = F16(x)
o=4*(x(1)^2)-2.1*(x(1)^4)+(x(1)^6)/3+x(1)*x(2)-4*(x(2)^2)+4*(x(2)^4);
end

% F17

function o = F17(x)
o=(x(2)-(x(1)^2)*5.1/(4*(pi^2))+5/pi*x(1)-6)^2+10*(1-1/(8*pi))*cos(x(1))+10;
end

% F18

function o = F18(x)
o=(1+(x(1)+x(2)+1)^2*(19-14*x(1)+3*(x(1)^2)-14*x(2)+6*x(1)*x(2)+3*x(2)^2))*...
    (30+(2*x(1)-3*x(2))^2*(18-32*x(1)+12*(x(1)^2)+48*x(2)-36*x(1)*x(2)+27*(x(2)^2)));
end

% F19

function o = F19(x)
aH=[3 10 30;.1 10 35;3 10 30;.1 10 35];cH=[1 1.2 3 3.2];
pH=[.3689 .117 .2673;.4699 .4387 .747;.1091 .8732 .5547;.03815 .5743 .8828];
o=0;
for i=1:4
    o=o-cH(i)*exp(-(sum(aH(i,:).*((x-pH(i,:)).^2))));
end
end

% F20

function o = F20(x)
aH=[10 3 17 3.5 1.7 8;.05 10 17 .1 8 14;3 3.5 1.7 10 17 8;17 8 .05 10 .1 14];
cH=[1 1.2 3 3.2];
pH=[.1312 .1696 .5569 .0124 .8283 .5886;.2329 .4135 .8307 .3736 .1004 .9991;...
.2348 .1415 .3522 .2883 .3047 .6650;.4047 .8828 .8732 .5743 .1091 .0381];
o=0;
for i=1:4
    o=o-cH(i)*exp(-(sum(aH(i,:).*((x-pH(i,:)).^2))));
end
end

% F21

function o = F21(x)
aSH=[4 4 4 4;1 1 1 1;8 8 8 8;6 6 6 6;3 7 3 7;2 9 2 9;5 5 3 3;8 1 8 1;6 2 6 2;7 3.6 7 3.6];
cSH=[.1 .2 .2 .4 .4 .6 .3 .7 .5 .5];

o=0;
for i=1:5
    o=o-((x-aSH(i,:))*(x-aSH(i,:))'+cSH(i))^(-1);
end
end

% F22

function o = F22(x)
aSH=[4 4 4 4;1 1 1 1;8 8 8 8;6 6 6 6;3 7 3 7;2 9 2 9;5 5 3 3;8 1 8 1;6 2 6 2;7 3.6 7 3.6];
cSH=[.1 .2 .2 .4 .4 .6 .3 .7 .5 .5];

o=0;
for i=1:7
    o=o-((x-aSH(i,:))*(x-aSH(i,:))'+cSH(i))^(-1);
end
end

% F23

function o = F23(x)
aSH=[4 4 4 4;1 1 1 1;8 8 8 8;6 6 6 6;3 7 3 7;2 9 2 9;5 5 3 3;8 1 8 1;6 2 6 2;7 3.6 7 3.6];
cSH=[.1 .2 .2 .4 .4 .6 .3 .7 .5 .5];

o=0;
for i=1:10
    o=o-((x-aSH(i,:))*(x-aSH(i,:))'+cSH(i))^(-1);
end
end

function o=Ufun(x,a,k,m)
o=k.*((x-a).^m).*(x>a)+k.*((-x-a).^m).*(x<(-a));
end
%% ====================CEC2014============================
function o = F24(x)
    o=cec14_func(x',1);
end
function o = F25(x)
    o=cec14_func(x',2);
end
function o = F26(x)
    o=cec14_func(x',3);
end
function o = F27(x)
    o=cec14_func(x',4);
end
function o = F28(x)
    o=cec14_func(x',5);
end
function o = F29(x)
    o=cec14_func(x',6);
end
function o = F30(x)
    o=cec14_func(x',7);
end
function o = F31(x)
    o=cec14_func(x',8);
end
function o = F32(x)
    o=cec14_func(x',9);
end
function o = F33(x)
    o=cec14_func(x',10);
end
function o = F34(x)
    o=cec14_func(x',11);
end
function o = F35(x)
    o=cec14_func(x',12);
end
function o = F36(x)
    o=cec14_func(x',13);
end
function o = F37(x)
    o=cec14_func(x',14);
end
function o = F38(x)
    o=cec14_func(x',15);
end
function o = F39(x)
    o=cec14_func(x',16);
end
function o = F40(x)
    o=cec14_func(x',17);
end
function o = F41(x)
    o=cec14_func(x',18);
end
function o = F42(x)
    o=cec14_func(x',19);
end
function o = F43(x)
    o=cec14_func(x',20);
end
function o = F44(x)
    o=cec14_func(x',21);
end
function o = F45(x)
    o=cec14_func(x',22);
end
function o = F46(x)
    o=cec14_func(x',23);
end
function o = F47(x)
    o=cec14_func(x',24);
end
function o = F48(x)
    o=cec14_func(x',25);
end
function o = F49(x)
    o=cec14_func(x',26);
end
function o = F50(x)
    o=cec14_func(x',27);
end
function o = F51(x)
    o=cec14_func(x',28);
end
function o = F52(x)
    o=cec14_func(x',29);
end
function o = F53(x)
    o=cec14_func(x',30);
end
%% ====================CEC2005============================
function o = F54(x)
    o=cec05_func(x,1);
end
function o = F55(x)
    o=cec05_func(x,2);
end
function o = F56(x)
    o=cec05_func(x,3);
end
function o = F57(x)
    o=cec05_func(x,4);
end
function o = F58(x)
    o=cec05_func(x,5);
end
function o = F59(x)
    o=cec05_func(x,6);
end
function o = F60(x)
    o=cec05_func(x,7);
end
function o = F61(x)
    o=cec05_func(x,8);
end
function o = F62(x)
    o=cec05_func(x,9);
end
function o = F63(x)
    o=cec05_func(x,10);
end
function o = F64(x)
    o=cec05_func(x,11);
end
function o = F65(x)
    o=cec05_func(x,12);
end
function o = F66(x)
    o=cec05_func(x,13);
end
function o = F67(x)
    o=cec05_func(x,14);
end
function o = F68(x)
    o=cec05_func(x,15);
end
function o = F69(x)
    o=cec05_func(x,16);
end
function o = F70(x)
    o=cec05_func(x,17);
end
function o = F71(x)
    o=cec05_func(x,18);
end
function o = F72(x)
    o=cec05_func(x,19);
end
function o = F73(x)
    o=cec05_func(x,20);
end
function o = F74(x)
    o=cec05_func(x,21);
end
function o = F75(x)
    o=cec05_func(x,22);
end
function o = F76(x)
    o=cec05_func(x,23);
end
function o = F77(x)
    o=cec05_func(x,24);
end
function o = F78(x)
    o=cec05_func(x,25);
end
%% ====================CEC2013============================
function o = F79(x)
    fhd=str2func('cec13_func');
    o=feval(fhd,x',1);
end
function o = F80(x)
    fhd=str2func('cec13_func');
    o=feval(fhd,x',2);
end
function o = F81(x)
  fhd=str2func('cec13_func');
    o=feval(fhd,x',3);
end
function o = F82(x)
   fhd=str2func('cec13_func');
    o=feval(fhd,x',4);
end
function o = F83(x)
    fhd=str2func('cec13_func');
    o=feval(fhd,x',5);
end
function o = F84(x)
    fhd=str2func('cec13_func');
    o=feval(fhd,x',6);
end
function o = F85(x)
    fhd=str2func('cec13_func');
    o=feval(fhd,x',7);
end
function o = F86(x)
    fhd=str2func('cec13_func');
    o=feval(fhd,x',8);
end
function o = F87(x)
    fhd=str2func('cec13_func');
    o=feval(fhd,x',9);
end
function o = F88(x)
    fhd=str2func('cec13_func');
    o=feval(fhd,x',10);
end
function o = F89(x)
    fhd=str2func('cec13_func');
    o=feval(fhd,x',11);
end
function o = F90(x)
   fhd=str2func('cec13_func');
    o=feval(fhd,x',12);
end
function o = F91(x)
    fhd=str2func('cec13_func');
    o=feval(fhd,x',13);
end
function o = F92(x)
    fhd=str2func('cec13_func');
    o=feval(fhd,x',14);
end
function o = F93(x)
    fhd=str2func('cec13_func');
    o=feval(fhd,x',15);
end
function o = F94(x)
    fhd=str2func('cec13_func');
    o=feval(fhd,x',16);
end
function o = F95(x)
    fhd=str2func('cec13_func');
    o=feval(fhd,x',17);
end
function o = F96(x)
    fhd=str2func('cec13_func');
    o=feval(fhd,x',18);
end
function o = F97(x)
    fhd=str2func('cec13_func');
    o=feval(fhd,x',19);
end
function o = F98(x)
    fhd=str2func('cec13_func');
    o=feval(fhd,x',20);
end
function o = F99(x)
    fhd=str2func('cec13_func');
    o=feval(fhd,x',21);
end
function o = F100(x)
    fhd=str2func('cec13_func');
    o=feval(fhd,x',22);
end
function o = F101(x)
    fhd=str2func('cec13_func');
    o=feval(fhd,x',23);
end
function o = F102(x)
    fhd=str2func('cec13_func');
    o=feval(fhd,x',24);
end
function o = F103(x)
    fhd=str2func('cec13_func');
    o=feval(fhd,x',25);
end
function o = F104(x)
    fhd=str2func('cec13_func');
    o=feval(fhd,x',26);
end
function o = F105(x)
    fhd=str2func('cec13_func');
    o=feval(fhd,x',27);
end
function o = F106(x)
    fhd=str2func('cec13_func');
    o=feval(fhd,x',28);
end
%% ====================CEC2017============================
function o = F107(x)
    fhd=str2func('cec17_func');
    o=feval(fhd,x',1);
end
function o = F108(x)
    fhd=str2func('cec17_func');
    o=feval(fhd,x',2);
end
function o = F109(x)
    fhd=str2func('cec17_func');
    o=feval(fhd,x',3);
end
function o = F110(x)
    fhd=str2func('cec17_func');
    o=feval(fhd,x',4);
end
function o = F111(x)
    fhd=str2func('cec17_func');
    o=feval(fhd,x',5);
end
function o = F112(x)
    fhd=str2func('cec17_func');
    o=feval(fhd,x',6);
end
function o = F113(x)
    fhd=str2func('cec17_func');
    o=feval(fhd,x',7);
end
function o = F114(x)
    fhd=str2func('cec17_func');
    o=feval(fhd,x',8);
end
function o = F115(x)
    fhd=str2func('cec17_func');
    o=feval(fhd,x',9);
end
function o = F116(x)
    fhd=str2func('cec17_func');
    o=feval(fhd,x',10);
end
function o = F117(x)
    fhd=str2func('cec17_func');
    o=feval(fhd,x',11);
end
function o = F118(x)
    fhd=str2func('cec17_func');
    o=feval(fhd,x',12);
end
function o = F119(x)
    fhd=str2func('cec17_func');
    o=feval(fhd,x',13);
end
function o = F120(x)
    fhd=str2func('cec17_func');
    o=feval(fhd,x',14);
end
function o = F121(x)
    fhd=str2func('cec17_func');
    o=feval(fhd,x',15);
end
function o = F122(x)
    fhd=str2func('cec17_func');
    o=feval(fhd,x',16);
end
function o = F123(x)
    fhd=str2func('cec17_func');
    o=feval(fhd,x',17);
end
function o = F124(x)
    fhd=str2func('cec17_func');
    o=feval(fhd,x',18);
end
function o = F125(x)
    fhd=str2func('cec17_func');
    o=feval(fhd,x',19);
end
function o = F126(x)
    fhd=str2func('cec17_func');
    o=feval(fhd,x',20);
end
function o = F127(x)
    fhd=str2func('cec17_func');
    o=feval(fhd,x',21);
end
function o = F128(x)
    fhd=str2func('cec17_func');
    o=feval(fhd,x',22);
end
function o = F129(x)
    fhd=str2func('cec17_func');
    o=feval(fhd,x',23);
end
function o = F130(x)
    fhd=str2func('cec17_func');
    o=feval(fhd,x',24);
end
function o = F131(x)
    fhd=str2func('cec17_func');
    o=feval(fhd,x',25);
end
function o = F132(x)
    fhd=str2func('cec17_func');
    o=feval(fhd,x',26);
end
function o = F133(x)
    fhd=str2func('cec17_func');
    o=feval(fhd,x',27);
end
function o = F134(x)
    fhd=str2func('cec17_func');
    o=feval(fhd,x',28);
end
function o = F135(x)
    fhd=str2func('cec17_func');
    o=feval(fhd,x',29);
end
function o = F136(x)
    fhd=str2func('cec17_func');
    o=feval(fhd,x',30);
end
%% ====================CEC2019============================
function o = F137(x)
    fhd=str2func('cec19_func');
    o=feval(fhd,x',1);
end
function o = F138(x)
    fhd=str2func('cec19_func');
    o=feval(fhd,x',2);
end
function o = F139(x)
    fhd=str2func('cec19_func');
    o=feval(fhd,x',3);
end
function o = F140(x)
    fhd=str2func('cec19_func');
    o=feval(fhd,x',4);
end
function o = F141(x)
    fhd=str2func('cec19_func');
    o=feval(fhd,x',5);
end
function o = F142(x)
    fhd=str2func('cec19_func');
    o=feval(fhd,x',6);
end
function o = F143(x)
    fhd=str2func('cec19_func');
    o=feval(fhd,x',7);
end
function o = F144(x)
    fhd=str2func('cec19_func');
    o=feval(fhd,x',8);
end
function o = F145(x)
    fhd=str2func('cec19_func');
    o=feval(fhd,x',9);
end
function o = F146(x)
    fhd=str2func('cec19_func');
    o=feval(fhd,x',10);
end