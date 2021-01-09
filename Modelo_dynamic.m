function [residual, g1, g2, g3] = Modelo_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Inputs :
%   y         [#dynamic variables by 1] double    vector of endogenous variables in the order stored
%                                                 in M_.lead_lag_incidence; see the Manual
%   x         [nperiods by M_.exo_nbr] double     matrix of exogenous variables (in declaration order)
%                                                 for all simulation periods
%   steady_state  [M_.endo_nbr by 1] double       vector of steady state values
%   params    [M_.param_nbr by 1] double          vector of parameter values in declaration order
%   it_       scalar double                       time period for exogenous variables for which to evaluate the model
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the dynamic model equations in order of 
%                                          declaration of the equations.
%                                          Dynare may prepend auxiliary equations, see M_.aux_vars
%   g1        [M_.endo_nbr by #dynamic variables] double    Jacobian matrix of the dynamic model equations;
%                                                           rows: equations in order of declaration
%                                                           columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g2        [M_.endo_nbr by (#dynamic variables)^2] double   Hessian matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g3        [M_.endo_nbr by (#dynamic variables)^3] double   Third order derivative matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(82, 1);
T9 = exp(y(42));
T11 = T9^params(62);
T16 = (1-exp(y(43)))^params(63);
T29 = (1-exp(y(73)))*exp(y(44))*exp(y(108))/exp(y(116));
T33 = exp(y(124))/T9;
T47 = T9^(params(62)/params(64));
T53 = exp(y(108))*(exp(y(125))-1)/(exp(y(125))*exp(y(115)));
T55 = T53^((-1)/params(64));
T62 = 1/exp(y(114));
T65 = exp(y(108))*T9^(-params(62));
T77 = T65*exp(y(107))*exp(y(103))*(1-exp(y(109))*params(83));
T82 = params(2)*exp(y(135))*exp(y(124))^(-params(62));
T103 = T62*(T77-T82*(exp(y(137))*exp(y(132))+exp(y(134))*exp(y(133))*params(75)*(1-params(83)*exp(y(136)))));
T118 = params(3)*exp(y(81));
T119 = exp(y(1))*exp(y(45))/T118;
T141 = exp(y(112))*exp(y(102))*exp(y(6))/params(3);
T169 = exp(y(101))*exp(y(103))*exp(y(107))*exp(y(109))*params(83);
T173 = exp(y(33))*params(84)/T118;
T183 = exp(y(33))*exp(y(110))/T118;
T227 = exp(y(127))*exp(y(126))*(exp(y(48))+params(5))/exp(y(130));
T238 = exp(y(7))*exp(y(90))/T118;
T241 = exp(y(4))*exp(y(54))/params(3);
T242 = (1-exp(y(74)))*T241;
T272 = params(61)*exp(y(131))/exp(y(125));
T314 = exp(y(27))*exp(y(93))/T118;
T333 = exp(y(96))/(1-exp(y(113)));
T334 = (-exp(y(113)))/(1-exp(y(113)))+T333;
T371 = exp(y(33))*exp(y(110))+exp(y(7))*exp(y(90));
T378 = exp(y(19))*params(51)+exp(y(27))*exp(y(93))+exp(y(1))*exp(y(45))+exp(y(120))*T371-exp(y(24))-params(87);
TEF_0 = logncdf(exp(y(119)),params(88),exp(y(121)));
T384 = TEF_0;
T394 = exp(y(47))*exp(y(50));
T400 = params(10)*exp(y(4))^(params(10)-1);
T402 = exp(y(5))^params(11);
T405 = exp(y(6))^params(12);
T410 = params(3)^(1-params(10)-params(11)-params(12));
T411 = T400*T402*T405*T410;
T415 = (exp(y(43))*exp(y(77)))^(1-params(10)-params(11)-params(12));
T423 = params(11)*exp(y(5))^(params(11)-1);
T424 = exp(y(4))^params(10);
T427 = T410*T405*T423*T424;
T435 = params(12)*exp(y(6))^(params(12)-1);
T438 = T410*T402*T424*T435;
T448 = params(3)^((-params(10))-params(11)-params(12));
T450 = exp(y(77))^(1-params(10)-params(11)-params(12));
T451 = T405*T402*(1-params(10)-params(11)-params(12))*T424*T448*T450;
T452 = exp(y(43))^((-params(10))-params(11)-params(12));
T459 = (exp(y(4))/params(3))^params(10);
T461 = (exp(y(5))/params(3))^params(11);
T464 = (exp(y(6))/params(3))^params(12);
T465 = T459*T461*T464;
T503 = (1-params(32))*exp(y(87))^(1-params(49))+params(32)*exp(y(8))^(1-params(49));
T516 = (params(80)-1)/params(80);
T524 = params(81)*exp(y(104))^T516+(1-params(81))*exp(y(105))^T516;
T534 = exp(y(103))*(1-params(81))/(params(81)*exp(y(106)));
T555 = (params(15)-1)/params(15);
T563 = params(16)*exp(y(62))^T555+(1-params(16))*exp(y(61))^T555;
T575 = (1-params(16))*exp(y(65))/(params(16)*exp(y(66)));
T602 = (params(24)-1)/params(24);
T610 = params(22)*exp(y(64))^T602+(1-params(22))*exp(y(68))^T602;
T622 = (1-params(22))*exp(y(70))/(params(22)*exp(y(100)));
T630 = (params(27)-1)/params(27);
T636 = params(25)*exp(y(64))^T630+(1-params(25))*exp(y(61))^T630;
T646 = (1-params(25))*exp(y(67))/(exp(y(66))*params(25));
T665 = exp(y(47))*exp(y(9))*(1+exp(y(72)))/params(3);
T673 = exp(y(6))*exp(y(102))*(exp(y(112))-1)/params(3);
T730 = exp(y(5))*exp(y(55))/(params(3)*exp(y(47)));
T735 = exp(y(71))/exp(y(47));
T746 = exp(y(9))*exp(y(72))/params(3);
T748 = exp(y(63))+exp(y(50))+exp(y(76))+exp(y(75))+exp(y(83))*exp(y(84))-exp(y(9))/params(3)-T746-T730;
T768 = exp(y(108))*T9^(1-params(62))/(1-params(62));
T772 = exp(y(116))*(1-exp(y(43)))^(1-params(63))/(1-params(63));
T777 = exp(y(114))*exp(y(53))^(1-params(65))/(1-params(65));
T782 = exp(y(115))*exp(y(80))^(1-params(64))/(1-params(64));
T796 = params(13)^(1-params(93));
T800 = exp(y(16))^params(93)*T796*exp(x(it_, 1));
T808 = params(8)^(1-params(94));
T812 = exp(y(14))^params(94)*T808*exp(x(it_, 2));
T820 = params(43)^(1-params(95));
T824 = exp(y(3))^params(95)*T820*exp(x(it_, 3));
T832 = params(46)^(1-params(96));
T833 = exp(y(21))^params(96)*T832;
T844 = params(47)^(1-params(97));
T848 = exp(y(22))^params(97)*T844*exp(x(it_, 5));
T856 = params(29)^(1-params(98));
T860 = exp(y(28))^params(98)*T856*exp(x(it_, 6));
T868 = params(9)^(1-params(99));
T872 = exp(y(15))^params(99)*T868*exp(x(it_, 7));
T880 = params(4)^(1-params(101));
T884 = exp(y(12))^params(101)*T880*exp(x(it_, 9));
T892 = params(17)^(1-params(102));
T896 = exp(y(17))^params(102)*T892*exp(x(it_, 10));
T904 = params(36)^(1-params(103));
T908 = exp(y(18))^params(103)*T904*exp(x(it_, 11));
T916 = params(7)^(1-params(104));
T920 = exp(y(13))^params(104)*T916*exp(x(it_, 12));
T928 = params(53)^(1-params(105));
T932 = exp(y(25))^params(105)*T928*exp(x(it_, 13));
T939 = params(42)^(1-params(106));
T943 = y(20)^params(106)*T939*exp(x(it_, 14));
T951 = params(28)^(1-params(107));
T955 = exp(y(10))^params(107)*T951*exp(x(it_, 15));
T963 = params(82)^(1-params(108));
T964 = exp(y(29))^params(108)*T963;
T975 = params(76)^(1-params(109));
T976 = exp(y(30))^params(109)*T975;
T987 = params(85)^(1-params(110));
T991 = exp(y(31))^params(110)*T987*exp(x(it_, 18));
T999 = params(77)^(1-params(111));
T1003 = exp(y(34))^params(111)*T999*exp(x(it_, 19));
T1011 = params(54)^(1-params(112));
T1015 = exp(y(35))^params(112)*T1011*exp(x(it_, 20));
T1023 = params(1)^(1-params(113));
T1027 = exp(y(38))^params(113)*T1023*exp(x(it_, 21));
T1035 = params(41)^(1-params(114));
T1039 = exp(y(37))^params(114)*T1035*exp(x(it_, 22));
T1047 = params(78)^(1-params(115));
T1048 = exp(y(36))^params(115)*T1047;
T1059 = params(90)^(1-params(116));
T1063 = exp(y(39))^params(116)*T1059*exp(x(it_, 24));
T1071 = params(89)^(1-params(117));
T1075 = exp(y(40))^params(117)*T1071*exp(x(it_, 25));
T1083 = params(86)^(1-params(118));
T1087 = exp(y(41))^params(118)*T1083*exp(x(it_, 26));
lhs =T11/T16;
rhs =T29;
residual(1)= lhs-rhs;
lhs =T33^params(62);
rhs =params(2)*exp(y(125))/exp(y(127));
residual(2)= lhs-rhs;
lhs =exp(y(80));
rhs =T47*T55;
residual(3)= lhs-rhs;
lhs =exp(y(53));
rhs =T103^((-1)/params(65));
residual(4)= lhs-rhs;
lhs =exp(y(46))-T119;
rhs =exp(y(111))-exp(y(33))/T118+(1-exp(y(73)))*(exp(y(43))*exp(y(44))+exp(y(86)))+T141+params(50)+exp(y(91))+exp(y(47))*(exp(y(75))+exp(y(76)))+exp(y(19))/T118-T9-exp(y(80))-exp(y(107))*exp(y(103))*exp(y(101))+T169-T173-exp(y(33))*(exp(y(110))-1)/T118;
residual(5)= lhs-rhs;
lhs =exp(y(111));
rhs =T169+T183-T173;
residual(6)= lhs-rhs;
lhs =exp(y(53));
rhs =exp(y(101))+params(75)*exp(y(6))/params(3);
residual(7)= lhs-rhs;
lhs =exp(y(48));
rhs =1/params(6)*(exp(y(127))*exp(y(126))/exp(y(130))-1);
residual(8)= lhs-rhs;
lhs =exp(y(48));
rhs =exp(y(49))/exp(y(4));
residual(9)= lhs-rhs;
lhs =exp(y(56));
rhs =(1-exp(y(74)))*exp(y(54))-params(6)/2*exp(y(48))^2-exp(y(48))+T227;
residual(10)= lhs-rhs;
lhs =exp(y(57));
rhs =params(50)+exp(y(49))+T238-T242-exp(y(47))*(1-params(48))*exp(y(83))*exp(y(84));
residual(11)= lhs-rhs;
lhs =exp(y(51));
rhs =exp(y(49))+exp(y(4))*params(5)/params(3);
residual(12)= lhs-rhs;
lhs =exp(y(123));
rhs =exp(y(111))+exp(y(57));
residual(13)= lhs-rhs;
lhs =exp(y(95))/exp(y(46));
rhs =T272^(-params(60));
residual(14)= lhs-rhs;
lhs =exp(y(26))*exp(y(96));
rhs =exp(y(1))*exp(y(45))+exp(y(27))*exp(y(93));
residual(15)= lhs-rhs;
lhs =exp(y(94));
rhs =exp(y(46))+exp(y(95));
residual(16)= lhs-rhs;
lhs =exp(y(92));
rhs =exp(y(113))*(exp(y(94))+exp(y(80))*params(51));
residual(17)= lhs-rhs;
lhs =exp(y(94))+exp(y(80))*params(51);
rhs =exp(y(111))+exp(y(57))+exp(y(92));
residual(18)= lhs-rhs;
lhs =exp(y(91));
rhs =T173+exp(y(80))*params(51)+exp(y(95))+exp(y(46))+T183+T238+exp(y(24))/T118-T119-T314-exp(y(57))-exp(y(111))-exp(y(92))-exp(y(19))*params(51)/T118-T169;
residual(19)= lhs-rhs;
lhs =exp(y(90));
rhs =params(91)*exp(y(23))+(1-params(91))*(T334+params(56)*params(55)*exp(y(57))^(params(55)-1)+exp(y(122))*exp(y(118)));
residual(20)= lhs-rhs;
lhs =exp(y(110));
rhs =params(92)*exp(y(32))+(1-params(92))*(exp(y(122))*exp(y(118))+T334+params(58)*params(57)*exp(y(111))^(params(57)-1));
residual(21)= lhs-rhs;
lhs =exp(y(119));
rhs =T378/T371;
residual(22)= lhs-rhs;
lhs =exp(y(118));
rhs =T384;
residual(23)= lhs-rhs;
lhs =exp(y(52));
rhs =params(5)*exp(y(5))/params(3)+T394;
residual(24)= lhs-rhs;
lhs =T411*T415;
rhs =exp(y(54))/exp(y(85));
residual(25)= lhs-rhs;
lhs =T415*T427;
rhs =exp(y(55))/exp(y(85));
residual(26)= lhs-rhs;
lhs =T415*T438;
rhs =exp(y(102))/exp(y(85));
residual(27)= lhs-rhs;
lhs =T451*T452;
rhs =exp(y(44))/exp(y(85));
residual(28)= lhs-rhs;
lhs =exp(y(59));
rhs =T415*T465;
residual(29)= lhs-rhs;
lhs =exp(y(87));
rhs =(1+params(31))*exp(y(88))/exp(y(89));
residual(30)= lhs-rhs;
lhs =exp(y(88));
rhs =exp(y(85))*exp(y(59))+params(2)*params(32)*exp(y(128));
residual(31)= lhs-rhs;
lhs =exp(y(89));
rhs =exp(y(59))+params(2)*params(32)*exp(y(129));
residual(32)= lhs-rhs;
lhs =exp(y(58));
rhs =T503^(1/(1-params(49)));
residual(33)= lhs-rhs;
lhs =exp(y(86));
rhs =exp(y(59))*(exp(y(58))-exp(y(85)));
residual(34)= lhs-rhs;
lhs =exp(y(59));
rhs =params(79)*T524^(params(80)/(params(80)-1));
residual(35)= lhs-rhs;
lhs =exp(y(104))/exp(y(105));
rhs =T534^(-params(80));
residual(36)= lhs-rhs;
lhs =exp(y(59))*exp(y(58));
rhs =exp(y(103))*exp(y(104))+exp(y(105))*exp(y(106));
residual(37)= lhs-rhs;
lhs =exp(y(104));
rhs =exp(y(101));
residual(38)= lhs-rhs;
lhs =T394+T9+exp(y(49))+exp(y(60));
rhs =params(14)*T563^(params(15)/(params(15)-1));
residual(39)= lhs-rhs;
lhs =exp(y(62))/exp(y(61));
rhs =T575^(-params(15));
residual(40)= lhs-rhs;
lhs =exp(y(65));
rhs =exp(y(47))*params(30)*(1+exp(y(79)));
residual(41)= lhs-rhs;
lhs =T394+T9+exp(y(49))+exp(y(60));
rhs =(1+exp(y(78)))*(exp(y(62))*exp(y(65))+exp(y(61))*exp(y(66)));
residual(42)= lhs-rhs;
lhs =exp(y(69));
rhs =params(23)*T610^(params(24)/(params(24)-1));
residual(43)= lhs-rhs;
lhs =exp(y(64))/exp(y(68));
rhs =T622^(-params(24));
residual(44)= lhs-rhs;
lhs =exp(y(105));
rhs =params(26)*T636^(params(27)/(params(27)-1));
residual(45)= lhs-rhs;
lhs =exp(y(64))/exp(y(61));
rhs =T646^(-params(27));
residual(46)= lhs-rhs;
lhs =exp(y(105))*exp(y(106));
rhs =exp(y(61))*exp(y(66))+exp(y(64))*exp(y(67));
residual(47)= lhs-rhs;
lhs =exp(y(67));
rhs =exp(y(47))*exp(y(70));
residual(48)= lhs-rhs;
lhs =exp(y(47))*exp(y(63));
rhs =exp(y(60))+T665-exp(y(74))*T241-exp(y(73))*(exp(y(43))*exp(y(44))+exp(y(86)))-T673-(T394+T9+exp(y(49))+exp(y(60)))*exp(y(78))/(1+exp(y(78)))-exp(y(62))*exp(y(47))*params(30)*exp(y(79))-exp(y(101))*exp(y(103))*(exp(y(107))-1)-exp(y(47))*exp(y(84))*params(48)*exp(y(83))-params(52);
residual(49)= lhs-rhs;
lhs =exp(y(72));
rhs =params(33)+exp(y(47))*exp(y(63))*params(34)/exp(y(71));
residual(50)= lhs-rhs;
lhs =exp(y(71));
rhs =exp(y(84))+T394+exp(y(64))+exp(y(49))+T9+exp(y(60))+exp(y(103))*exp(y(101))-exp(y(62))/(1+exp(y(79)));
residual(51)= lhs-rhs;
lhs =y(97);
rhs =exp(y(50))+exp(y(63))+exp(y(76))+exp(y(75))+exp(y(64))*exp(y(70))+exp(y(83))*exp(y(84))-y(82)-exp(y(62))*params(30)-exp(y(9))*(1+exp(y(72)))/params(3)-T730;
residual(52)= lhs-rhs;
lhs =y(98);
rhs =(exp(y(64))*exp(y(70))-exp(y(62))*params(30))/T735;
residual(53)= lhs-rhs;
lhs =y(99);
rhs =T748/T735;
residual(54)= lhs-rhs;
lhs =exp(y(80))*(1-params(51));
rhs =exp(y(95))+exp(y(24))/T118+params(52)+exp(y(19))*(1-params(51))/T118+exp(y(47))*y(82)-exp(y(92))-T314;
residual(55)= lhs-rhs;
lhs =exp(y(117));
rhs =T768+T772+T777+T782+params(2)*params(3)*exp(y(138));
residual(56)= lhs-rhs;
lhs =exp(y(77));
rhs =T800;
residual(57)= lhs-rhs;
lhs =exp(y(75));
rhs =T812;
residual(58)= lhs-rhs;
lhs =exp(y(50));
rhs =T824;
residual(59)= lhs-rhs;
lhs =exp(y(83));
rhs =T833/exp(x(it_, 4));
residual(60)= lhs-rhs;
lhs =exp(y(84));
rhs =T848;
residual(61)= lhs-rhs;
lhs =exp(y(100));
rhs =T860;
residual(62)= lhs-rhs;
lhs =exp(y(76));
rhs =T872;
residual(63)= lhs-rhs;
lhs =exp(y(73));
rhs =T884;
residual(64)= lhs-rhs;
lhs =exp(y(78));
rhs =T896;
residual(65)= lhs-rhs;
lhs =exp(y(79));
rhs =T908;
residual(66)= lhs-rhs;
lhs =exp(y(74));
rhs =T920;
residual(67)= lhs-rhs;
lhs =exp(y(93));
rhs =T932;
residual(68)= lhs-rhs;
lhs =y(82);
rhs =T943;
residual(69)= lhs-rhs;
lhs =exp(y(69));
rhs =T955;
residual(70)= lhs-rhs;
lhs =exp(y(107));
rhs =T964/exp(x(it_, 16));
residual(71)= lhs-rhs;
lhs =exp(y(108));
rhs =T976/exp(x(it_, 17));
residual(72)= lhs-rhs;
lhs =exp(y(109));
rhs =T991;
residual(73)= lhs-rhs;
lhs =exp(y(112));
rhs =T1003;
residual(74)= lhs-rhs;
lhs =exp(y(113));
rhs =T1015;
residual(75)= lhs-rhs;
lhs =exp(y(116));
rhs =T1027;
residual(76)= lhs-rhs;
lhs =exp(y(115));
rhs =T1039;
residual(77)= lhs-rhs;
lhs =exp(y(114));
rhs =T1048/exp(x(it_, 23));
residual(78)= lhs-rhs;
lhs =exp(y(120));
rhs =T1063;
residual(79)= lhs-rhs;
lhs =exp(y(121));
rhs =T1075;
residual(80)= lhs-rhs;
lhs =exp(y(122));
rhs =T1087;
residual(81)= lhs-rhs;
lhs =exp(y(60));
rhs =params(18)-params(45)*(exp(y(9))*exp(y(2))/exp(y(11))-params(35));
residual(82)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(82, 164);

  %
  % Jacobian matrix
  %

T1109 = getPowerDeriv(T33,params(62),1);
T1120 = getPowerDeriv(T103,(-1)/params(65),1);
T1151 = exp(y(43))*exp(y(77))*getPowerDeriv(exp(y(43))*exp(y(77)),1-params(10)-params(11)-params(12),1);
T1152 = T411*T1151;
T1172 = (-(exp(y(1))*exp(y(45))/T371));
T1180 = getPowerDeriv(T53,(-1)/params(64),1);
T1189 = getPowerDeriv(T272,(-params(60)),1);
T1214 = (-(exp(y(5))*exp(y(55))*params(3)*exp(y(47))))/(params(3)*exp(y(47))*params(3)*exp(y(47)));
T1219 = (-(exp(y(47))*exp(y(71))))/(exp(y(47))*exp(y(47)));
T1270 = exp(y(4))*getPowerDeriv(exp(y(4)),params(10),1);
T1294 = exp(y(5))*getPowerDeriv(exp(y(5)),params(11),1);
T1323 = (-T141);
T1326 = exp(y(6))*getPowerDeriv(exp(y(6)),params(12),1);
T1359 = (-(1/params(6)*exp(y(127))*exp(y(126))/exp(y(130))));
T1368 = (-((T371*exp(y(7))*exp(y(90))*exp(y(120))-exp(y(7))*exp(y(90))*T378)/(T371*T371)));
T1378 = getPowerDeriv(T503,1/(1-params(49)),1);
T1394 = getPowerDeriv(T563,params(15)/(params(15)-1),1);
T1407 = getPowerDeriv(T636,params(27)/(params(27)-1),1);
T1440 = getPowerDeriv(T610,params(24)/(params(24)-1),1);
T1455 = getPowerDeriv(T575,(-params(15)),1);
T1468 = getPowerDeriv(T646,(-params(27)),1);
T1488 = getPowerDeriv(T622,(-params(24)),1);
T1622 = (-(T118*exp(y(33))*params(84)))/(T118*T118);
T1785 = getPowerDeriv(T534,(-params(80)),1);
T1797 = getPowerDeriv(T524,params(80)/(params(80)-1),1);
T1861 = (-((T371*exp(y(33))*exp(y(110))*exp(y(120))-exp(y(33))*exp(y(110))*T378)/(T371*T371)));
T1894 = ((-exp(y(113)))*(1-exp(y(113)))-(-exp(y(113)))*(-exp(y(113))))/((1-exp(y(113)))*(1-exp(y(113))))+(-(exp(y(96))*(-exp(y(113)))))/((1-exp(y(113)))*(1-exp(y(113))));
TEFD_fdd_0_1 = jacob_element('logncdf',1,{exp(y(119)),params(88),exp(y(121))});
T1938 = TEFD_fdd_0_1;
TEFD_fdd_0_3 = jacob_element('logncdf',3,{exp(y(119)),params(88),exp(y(121))});
T1941 = TEFD_fdd_0_3;
  g1(1,42)=T9*getPowerDeriv(T9,params(62),1)/T16;
  g1(1,43)=(-(T11*(-exp(y(43)))*getPowerDeriv(1-exp(y(43)),params(63),1)))/(T16*T16);
  g1(1,44)=(-T29);
  g1(1,73)=(-(exp(y(108))*exp(y(44))*(-exp(y(73)))/exp(y(116))));
  g1(1,108)=(-T29);
  g1(1,116)=(-((-((1-exp(y(73)))*exp(y(44))*exp(y(108))*exp(y(116))))/(exp(y(116))*exp(y(116)))));
  g1(2,42)=(-(T9*exp(y(124))))/(T9*T9)*T1109;
  g1(2,124)=T33*T1109;
  g1(2,125)=(-(params(2)*exp(y(125))/exp(y(127))));
  g1(2,127)=(-((-(params(2)*exp(y(125))*exp(y(127))))/(exp(y(127))*exp(y(127)))));
  g1(3,42)=(-(T55*T9*getPowerDeriv(T9,params(62)/params(64),1)));
  g1(3,125)=(-(T47*(exp(y(125))*exp(y(115))*exp(y(108))*exp(y(125))-exp(y(108))*(exp(y(125))-1)*exp(y(125))*exp(y(115)))/(exp(y(125))*exp(y(115))*exp(y(125))*exp(y(115)))*T1180));
  g1(3,80)=exp(y(80));
  g1(3,108)=(-(T47*T53*T1180));
  g1(3,115)=(-(T47*T1180*(-(exp(y(108))*(exp(y(125))-1)*exp(y(125))*exp(y(115))))/(exp(y(125))*exp(y(115))*exp(y(125))*exp(y(115)))));
  g1(4,42)=(-(T62*exp(y(107))*exp(y(103))*(1-exp(y(109))*params(83))*exp(y(108))*T9*getPowerDeriv(T9,(-params(62)),1)*T1120));
  g1(4,124)=(-(T1120*T62*(-((exp(y(137))*exp(y(132))+exp(y(134))*exp(y(133))*params(75)*(1-params(83)*exp(y(136))))*params(2)*exp(y(135))*exp(y(124))*getPowerDeriv(exp(y(124)),(-params(62)),1)))));
  g1(4,53)=exp(y(53));
  g1(4,132)=(-(T1120*T62*(-(T82*exp(y(137))*exp(y(132))))));
  g1(4,103)=(-(T1120*T62*T77));
  g1(4,133)=(-(T1120*T62*(-(T82*exp(y(134))*exp(y(133))*params(75)*(1-params(83)*exp(y(136)))))));
  g1(4,107)=(-(T1120*T62*T77));
  g1(4,134)=(-(T1120*T62*(-(T82*exp(y(134))*exp(y(133))*params(75)*(1-params(83)*exp(y(136)))))));
  g1(4,108)=(-(T1120*T62*T77));
  g1(4,135)=(-(T1120*T62*(-(T82*(exp(y(137))*exp(y(132))+exp(y(134))*exp(y(133))*params(75)*(1-params(83)*exp(y(136))))))));
  g1(4,109)=(-(T1120*T62*T65*exp(y(107))*exp(y(103))*(-(exp(y(109))*params(83)))));
  g1(4,136)=(-(T1120*T62*(-(T82*exp(y(134))*exp(y(133))*params(75)*(-(params(83)*exp(y(136))))))));
  g1(4,137)=(-(T1120*T62*(-(T82*exp(y(137))*exp(y(132))))));
  g1(4,114)=(-(T1120*(T77-T82*(exp(y(137))*exp(y(132))+exp(y(134))*exp(y(133))*params(75)*(1-params(83)*exp(y(136)))))*(-exp(y(114)))/(exp(y(114))*exp(y(114)))));
  g1(5,42)=T9;
  g1(5,43)=(-((1-exp(y(73)))*exp(y(43))*exp(y(44))));
  g1(5,44)=(-((1-exp(y(73)))*exp(y(43))*exp(y(44))));
  g1(5,45)=(-T119);
  g1(5,1)=(-T119);
  g1(5,46)=exp(y(46));
  g1(5,47)=(-(exp(y(47))*(exp(y(75))+exp(y(76)))));
  g1(5,6)=T1323;
  g1(5,73)=(-((exp(y(43))*exp(y(44))+exp(y(86)))*(-exp(y(73)))));
  g1(5,75)=(-(exp(y(47))*exp(y(75))));
  g1(5,76)=(-(exp(y(47))*exp(y(76))));
  g1(5,19)=(-(exp(y(19))/T118));
  g1(5,80)=exp(y(80));
  g1(5,81)=(-((-(exp(y(1))*exp(y(45))*T118))/(T118*T118)))-((-((-(T118*exp(y(33))))/(T118*T118)))+(-(T118*exp(y(19))))/(T118*T118)-T1622-(-(T118*exp(y(33))*(exp(y(110))-1)))/(T118*T118));
  g1(5,86)=(-((1-exp(y(73)))*exp(y(86))));
  g1(5,91)=(-exp(y(91)));
  g1(5,101)=(-(T169+(-(exp(y(107))*exp(y(103))*exp(y(101))))));
  g1(5,102)=T1323;
  g1(5,103)=(-(T169+(-(exp(y(107))*exp(y(103))*exp(y(101))))));
  g1(5,107)=(-(T169+(-(exp(y(107))*exp(y(103))*exp(y(101))))));
  g1(5,109)=(-T169);
  g1(5,110)=T183;
  g1(5,33)=(-((-(exp(y(33))/T118))-T173-exp(y(33))*(exp(y(110))-1)/T118));
  g1(5,111)=(-exp(y(111)));
  g1(5,112)=T1323;
  g1(6,81)=(-((-(T118*exp(y(33))*exp(y(110))))/(T118*T118)-T1622));
  g1(6,101)=(-T169);
  g1(6,103)=(-T169);
  g1(6,107)=(-T169);
  g1(6,109)=(-T169);
  g1(6,110)=(-T183);
  g1(6,33)=(-(T183-T173));
  g1(6,111)=exp(y(111));
  g1(7,6)=(-(params(75)*exp(y(6))/params(3)));
  g1(7,53)=exp(y(53));
  g1(7,101)=(-exp(y(101)));
  g1(8,48)=exp(y(48));
  g1(8,126)=T1359;
  g1(8,127)=T1359;
  g1(8,130)=(-(1/params(6)*(-(exp(y(127))*exp(y(126))*exp(y(130))))/(exp(y(130))*exp(y(130)))));
  g1(9,48)=exp(y(48));
  g1(9,49)=(-(exp(y(49))/exp(y(4))));
  g1(9,4)=(-((-(exp(y(49))*exp(y(4))))/(exp(y(4))*exp(y(4)))));
  g1(10,48)=(-((-(params(6)/2*exp(y(48))*2*exp(y(48))))-exp(y(48))+exp(y(127))*exp(y(48))*exp(y(126))/exp(y(130))));
  g1(10,54)=(-((1-exp(y(74)))*exp(y(54))));
  g1(10,56)=exp(y(56));
  g1(10,126)=(-T227);
  g1(10,74)=(-(exp(y(54))*(-exp(y(74)))));
  g1(10,127)=(-T227);
  g1(10,130)=(-((-(exp(y(130))*exp(y(127))*exp(y(126))*(exp(y(48))+params(5))))/(exp(y(130))*exp(y(130)))));
  g1(11,47)=exp(y(47))*(1-params(48))*exp(y(83))*exp(y(84));
  g1(11,49)=(-exp(y(49)));
  g1(11,4)=T242;
  g1(11,54)=T242;
  g1(11,7)=(-T238);
  g1(11,57)=exp(y(57));
  g1(11,74)=T241*(-exp(y(74)));
  g1(11,81)=(-((-(T118*exp(y(7))*exp(y(90))))/(T118*T118)));
  g1(11,83)=exp(y(47))*(1-params(48))*exp(y(83))*exp(y(84));
  g1(11,84)=exp(y(47))*(1-params(48))*exp(y(83))*exp(y(84));
  g1(11,90)=(-T238);
  g1(12,49)=(-exp(y(49)));
  g1(12,4)=(-(exp(y(4))*params(5)/params(3)));
  g1(12,51)=exp(y(51));
  g1(13,57)=(-exp(y(57)));
  g1(13,111)=(-exp(y(111)));
  g1(13,123)=exp(y(123));
  g1(14,125)=(-(params(61)*(-(exp(y(125))*exp(y(131))))/(exp(y(125))*exp(y(125)))*T1189));
  g1(14,46)=(-(exp(y(46))*exp(y(95))))/(exp(y(46))*exp(y(46)));
  g1(14,131)=(-(T272*T1189));
  g1(14,95)=exp(y(95))/exp(y(46));
  g1(15,45)=(-(exp(y(1))*exp(y(45))));
  g1(15,1)=(-(exp(y(1))*exp(y(45))));
  g1(15,93)=(-(exp(y(27))*exp(y(93))));
  g1(15,26)=exp(y(26))*exp(y(96));
  g1(15,27)=(-(exp(y(27))*exp(y(93))));
  g1(15,96)=exp(y(26))*exp(y(96));
  g1(16,46)=(-exp(y(46)));
  g1(16,94)=exp(y(94));
  g1(16,95)=(-exp(y(95)));
  g1(17,80)=(-(exp(y(113))*exp(y(80))*params(51)));
  g1(17,92)=exp(y(92));
  g1(17,94)=(-(exp(y(94))*exp(y(113))));
  g1(17,113)=(-(exp(y(113))*(exp(y(94))+exp(y(80))*params(51))));
  g1(18,57)=(-exp(y(57)));
  g1(18,80)=exp(y(80))*params(51);
  g1(18,92)=(-exp(y(92)));
  g1(18,94)=exp(y(94));
  g1(18,111)=(-exp(y(111)));
  g1(19,45)=T119;
  g1(19,1)=T119;
  g1(19,46)=(-exp(y(46)));
  g1(19,7)=(-T238);
  g1(19,57)=exp(y(57));
  g1(19,19)=exp(y(19))*params(51)/T118;
  g1(19,80)=(-(exp(y(80))*params(51)));
  g1(19,81)=(-(T1622+(-(T118*exp(y(33))*exp(y(110))))/(T118*T118)+(-(T118*exp(y(7))*exp(y(90))))/(T118*T118)+(-(T118*exp(y(24))))/(T118*T118)-(-(exp(y(1))*exp(y(45))*T118))/(T118*T118)-(-(T118*exp(y(27))*exp(y(93))))/(T118*T118)-(-(T118*exp(y(19))*params(51)))/(T118*T118)));
  g1(19,90)=(-T238);
  g1(19,91)=exp(y(91));
  g1(19,24)=(-(exp(y(24))/T118));
  g1(19,92)=exp(y(92));
  g1(19,93)=T314;
  g1(19,27)=T314;
  g1(19,95)=(-exp(y(95)));
  g1(19,101)=T169;
  g1(19,103)=T169;
  g1(19,107)=T169;
  g1(19,109)=T169;
  g1(19,110)=(-T183);
  g1(19,33)=(-(T173+T183));
  g1(19,111)=exp(y(111));
  g1(20,57)=(-((1-params(91))*params(56)*params(55)*exp(y(57))*getPowerDeriv(exp(y(57)),params(55)-1,1)));
  g1(20,23)=(-(params(91)*exp(y(23))));
  g1(20,90)=exp(y(90));
  g1(20,96)=(-((1-params(91))*T333));
  g1(20,113)=(-((1-params(91))*T1894));
  g1(20,118)=(-((1-params(91))*exp(y(122))*exp(y(118))));
  g1(20,122)=(-((1-params(91))*exp(y(122))*exp(y(118))));
  g1(21,96)=(-(T333*(1-params(92))));
  g1(21,32)=(-(params(92)*exp(y(32))));
  g1(21,110)=exp(y(110));
  g1(21,111)=(-((1-params(92))*params(58)*params(57)*exp(y(111))*getPowerDeriv(exp(y(111)),params(57)-1,1)));
  g1(21,113)=(-((1-params(92))*T1894));
  g1(21,118)=(-(exp(y(122))*exp(y(118))*(1-params(92))));
  g1(21,122)=(-(exp(y(122))*exp(y(118))*(1-params(92))));
  g1(22,45)=T1172;
  g1(22,1)=T1172;
  g1(22,7)=T1368;
  g1(22,19)=(-(exp(y(19))*params(51)/T371));
  g1(22,90)=T1368;
  g1(22,24)=(-((-exp(y(24)))/T371));
  g1(22,93)=(-(exp(y(27))*exp(y(93))/T371));
  g1(22,27)=(-(exp(y(27))*exp(y(93))/T371));
  g1(22,110)=T1861;
  g1(22,33)=T1861;
  g1(22,119)=exp(y(119));
  g1(22,120)=(-(exp(y(120))*T371/T371));
  g1(23,118)=exp(y(118));
  g1(23,119)=(-(exp(y(119))*T1938));
  g1(23,121)=(-(exp(y(121))*T1941));
  g1(24,47)=(-T394);
  g1(24,50)=(-T394);
  g1(24,5)=(-(params(5)*exp(y(5))/params(3)));
  g1(24,52)=exp(y(52));
  g1(25,43)=T1152;
  g1(25,4)=T415*T410*T405*T402*params(10)*exp(y(4))*getPowerDeriv(exp(y(4)),params(10)-1,1);
  g1(25,5)=T415*T410*T405*T400*T1294;
  g1(25,6)=T415*T410*T400*T402*T1326;
  g1(25,54)=(-(exp(y(54))/exp(y(85))));
  g1(25,77)=T1152;
  g1(25,85)=(-((-(exp(y(54))*exp(y(85))))/(exp(y(85))*exp(y(85)))));
  g1(26,43)=T427*T1151;
  g1(26,4)=T415*T410*T405*T423*T1270;
  g1(26,5)=T415*T410*T405*T424*params(11)*exp(y(5))*getPowerDeriv(exp(y(5)),params(11)-1,1);
  g1(26,6)=T415*T410*T423*T424*T1326;
  g1(26,55)=(-(exp(y(55))/exp(y(85))));
  g1(26,77)=T427*T1151;
  g1(26,85)=(-((-(exp(y(85))*exp(y(55))))/(exp(y(85))*exp(y(85)))));
  g1(27,43)=T438*T1151;
  g1(27,4)=T415*T410*T402*T435*T1270;
  g1(27,5)=T415*T410*T424*T435*T1294;
  g1(27,6)=T415*T410*T402*T424*params(12)*exp(y(6))*getPowerDeriv(exp(y(6)),params(12)-1,1);
  g1(27,77)=T438*T1151;
  g1(27,85)=(-((-(exp(y(102))*exp(y(85))))/(exp(y(85))*exp(y(85)))));
  g1(27,102)=(-(exp(y(102))/exp(y(85))));
  g1(28,43)=T451*exp(y(43))*getPowerDeriv(exp(y(43)),(-params(10))-params(11)-params(12),1);
  g1(28,44)=(-(exp(y(44))/exp(y(85))));
  g1(28,4)=T452*T450*T448*T405*T402*(1-params(10)-params(11)-params(12))*T1270;
  g1(28,5)=T452*T450*T448*T405*(1-params(10)-params(11)-params(12))*T424*T1294;
  g1(28,6)=T452*T450*T448*T402*(1-params(10)-params(11)-params(12))*T424*T1326;
  g1(28,77)=T452*T405*T402*(1-params(10)-params(11)-params(12))*T424*T448*exp(y(77))*getPowerDeriv(exp(y(77)),1-params(10)-params(11)-params(12),1);
  g1(28,85)=(-((-(exp(y(44))*exp(y(85))))/(exp(y(85))*exp(y(85)))));
  g1(29,43)=(-(T465*T1151));
  g1(29,4)=(-(T415*T464*T461*exp(y(4))/params(3)*getPowerDeriv(exp(y(4))/params(3),params(10),1)));
  g1(29,5)=(-(T415*T464*T459*exp(y(5))/params(3)*getPowerDeriv(exp(y(5))/params(3),params(11),1)));
  g1(29,6)=(-(T415*T459*T461*exp(y(6))/params(3)*getPowerDeriv(exp(y(6))/params(3),params(12),1)));
  g1(29,59)=exp(y(59));
  g1(29,77)=(-(T465*T1151));
  g1(30,87)=exp(y(87));
  g1(30,88)=(-((1+params(31))*exp(y(88))/exp(y(89))));
  g1(30,89)=(-((1+params(31))*(-(exp(y(88))*exp(y(89))))/(exp(y(89))*exp(y(89)))));
  g1(31,59)=(-(exp(y(85))*exp(y(59))));
  g1(31,85)=(-(exp(y(85))*exp(y(59))));
  g1(31,88)=exp(y(88));
  g1(31,128)=(-(params(2)*params(32)*exp(y(128))));
  g1(32,59)=(-exp(y(59)));
  g1(32,89)=exp(y(89));
  g1(32,129)=(-(params(2)*params(32)*exp(y(129))));
  g1(33,8)=(-(params(32)*exp(y(8))*getPowerDeriv(exp(y(8)),1-params(49),1)*T1378));
  g1(33,58)=exp(y(58));
  g1(33,87)=(-(T1378*(1-params(32))*exp(y(87))*getPowerDeriv(exp(y(87)),1-params(49),1)));
  g1(34,58)=(-(exp(y(59))*exp(y(58))));
  g1(34,59)=(-(exp(y(59))*(exp(y(58))-exp(y(85)))));
  g1(34,85)=(-(exp(y(59))*(-exp(y(85)))));
  g1(34,86)=exp(y(86));
  g1(35,59)=exp(y(59));
  g1(35,104)=(-(params(79)*params(81)*exp(y(104))*getPowerDeriv(exp(y(104)),T516,1)*T1797));
  g1(35,105)=(-(params(79)*T1797*(1-params(81))*exp(y(105))*getPowerDeriv(exp(y(105)),T516,1)));
  g1(36,103)=(-(T534*T1785));
  g1(36,104)=exp(y(104))/exp(y(105));
  g1(36,105)=(-(exp(y(104))*exp(y(105))))/(exp(y(105))*exp(y(105)));
  g1(36,106)=(-(T1785*(-(exp(y(103))*(1-params(81))*params(81)*exp(y(106))))/(params(81)*exp(y(106))*params(81)*exp(y(106)))));
  g1(37,58)=exp(y(59))*exp(y(58));
  g1(37,59)=exp(y(59))*exp(y(58));
  g1(37,103)=(-(exp(y(103))*exp(y(104))));
  g1(37,104)=(-(exp(y(103))*exp(y(104))));
  g1(37,105)=(-(exp(y(105))*exp(y(106))));
  g1(37,106)=(-(exp(y(105))*exp(y(106))));
  g1(38,101)=(-exp(y(101)));
  g1(38,104)=exp(y(104));
  g1(39,42)=T9;
  g1(39,47)=T394;
  g1(39,49)=exp(y(49));
  g1(39,50)=T394;
  g1(39,60)=exp(y(60));
  g1(39,61)=(-(params(14)*(1-params(16))*exp(y(61))*getPowerDeriv(exp(y(61)),T555,1)*T1394));
  g1(39,62)=(-(params(14)*T1394*params(16)*exp(y(62))*getPowerDeriv(exp(y(62)),T555,1)));
  g1(40,61)=(-(exp(y(62))*exp(y(61))))/(exp(y(61))*exp(y(61)));
  g1(40,62)=exp(y(62))/exp(y(61));
  g1(40,65)=(-(T575*T1455));
  g1(40,66)=(-(T1455*(-((1-params(16))*exp(y(65))*params(16)*exp(y(66))))/(params(16)*exp(y(66))*params(16)*exp(y(66)))));
  g1(41,47)=(-(exp(y(47))*params(30)*(1+exp(y(79)))));
  g1(41,65)=exp(y(65));
  g1(41,79)=(-(exp(y(47))*params(30)*exp(y(79))));
  g1(42,42)=T9;
  g1(42,47)=T394;
  g1(42,49)=exp(y(49));
  g1(42,50)=T394;
  g1(42,60)=exp(y(60));
  g1(42,61)=(-((1+exp(y(78)))*exp(y(61))*exp(y(66))));
  g1(42,62)=(-((1+exp(y(78)))*exp(y(62))*exp(y(65))));
  g1(42,65)=(-((1+exp(y(78)))*exp(y(62))*exp(y(65))));
  g1(42,66)=(-((1+exp(y(78)))*exp(y(61))*exp(y(66))));
  g1(42,78)=(-(exp(y(78))*(exp(y(62))*exp(y(65))+exp(y(61))*exp(y(66)))));
  g1(43,64)=(-(params(23)*params(22)*exp(y(64))*getPowerDeriv(exp(y(64)),T602,1)*T1440));
  g1(43,68)=(-(params(23)*T1440*(1-params(22))*exp(y(68))*getPowerDeriv(exp(y(68)),T602,1)));
  g1(43,69)=exp(y(69));
  g1(44,64)=exp(y(64))/exp(y(68));
  g1(44,68)=(-(exp(y(64))*exp(y(68))))/(exp(y(68))*exp(y(68)));
  g1(44,70)=(-(T622*T1488));
  g1(44,100)=(-(T1488*(-((1-params(22))*exp(y(70))*params(22)*exp(y(100))))/(params(22)*exp(y(100))*params(22)*exp(y(100)))));
  g1(45,61)=(-(params(26)*(1-params(25))*exp(y(61))*getPowerDeriv(exp(y(61)),T630,1)*T1407));
  g1(45,64)=(-(params(26)*T1407*params(25)*exp(y(64))*getPowerDeriv(exp(y(64)),T630,1)));
  g1(45,105)=exp(y(105));
  g1(46,61)=(-(exp(y(61))*exp(y(64))))/(exp(y(61))*exp(y(61)));
  g1(46,64)=exp(y(64))/exp(y(61));
  g1(46,66)=(-((-((1-params(25))*exp(y(67))*exp(y(66))*params(25)))/(exp(y(66))*params(25)*exp(y(66))*params(25))*T1468));
  g1(46,67)=(-(T646*T1468));
  g1(47,61)=(-(exp(y(61))*exp(y(66))));
  g1(47,64)=(-(exp(y(64))*exp(y(67))));
  g1(47,66)=(-(exp(y(61))*exp(y(66))));
  g1(47,67)=(-(exp(y(64))*exp(y(67))));
  g1(47,105)=exp(y(105))*exp(y(106));
  g1(47,106)=exp(y(105))*exp(y(106));
  g1(48,47)=(-(exp(y(47))*exp(y(70))));
  g1(48,67)=exp(y(67));
  g1(48,70)=(-(exp(y(47))*exp(y(70))));
  g1(49,42)=T9*exp(y(78))/(1+exp(y(78)));
  g1(49,43)=exp(y(73))*exp(y(43))*exp(y(44));
  g1(49,44)=exp(y(73))*exp(y(43))*exp(y(44));
  g1(49,47)=exp(y(47))*exp(y(63))-(T665-T394*exp(y(78))/(1+exp(y(78)))-exp(y(62))*exp(y(47))*params(30)*exp(y(79))-exp(y(47))*exp(y(84))*params(48)*exp(y(83)));
  g1(49,49)=exp(y(49))*exp(y(78))/(1+exp(y(78)));
  g1(49,50)=T394*exp(y(78))/(1+exp(y(78)));
  g1(49,4)=exp(y(74))*T241;
  g1(49,6)=T673;
  g1(49,54)=exp(y(74))*T241;
  g1(49,60)=(-(exp(y(60))-exp(y(60))*exp(y(78))/(1+exp(y(78)))));
  g1(49,62)=exp(y(62))*exp(y(47))*params(30)*exp(y(79));
  g1(49,9)=(-T665);
  g1(49,63)=exp(y(47))*exp(y(63));
  g1(49,72)=(-(exp(y(47))*exp(y(9))*exp(y(72))/params(3)));
  g1(49,73)=exp(y(73))*(exp(y(43))*exp(y(44))+exp(y(86)));
  g1(49,74)=exp(y(74))*T241;
  g1(49,78)=((1+exp(y(78)))*(T394+T9+exp(y(49))+exp(y(60)))*exp(y(78))-exp(y(78))*(T394+T9+exp(y(49))+exp(y(60)))*exp(y(78)))/((1+exp(y(78)))*(1+exp(y(78))));
  g1(49,79)=exp(y(62))*exp(y(47))*params(30)*exp(y(79));
  g1(49,83)=exp(y(47))*exp(y(84))*params(48)*exp(y(83));
  g1(49,84)=exp(y(47))*exp(y(84))*params(48)*exp(y(83));
  g1(49,86)=exp(y(73))*exp(y(86));
  g1(49,101)=exp(y(101))*exp(y(103))*(exp(y(107))-1);
  g1(49,102)=T673;
  g1(49,103)=exp(y(101))*exp(y(103))*(exp(y(107))-1);
  g1(49,107)=exp(y(107))*exp(y(103))*exp(y(101));
  g1(49,112)=T141;
  g1(50,47)=(-(exp(y(47))*exp(y(63))*params(34)/exp(y(71))));
  g1(50,63)=(-(exp(y(47))*exp(y(63))*params(34)/exp(y(71))));
  g1(50,71)=(-((-(exp(y(47))*exp(y(63))*params(34)*exp(y(71))))/(exp(y(71))*exp(y(71)))));
  g1(50,72)=exp(y(72));
  g1(51,42)=(-T9);
  g1(51,47)=(-T394);
  g1(51,49)=(-exp(y(49)));
  g1(51,50)=(-T394);
  g1(51,60)=(-exp(y(60)));
  g1(51,62)=exp(y(62))/(1+exp(y(79)));
  g1(51,64)=(-exp(y(64)));
  g1(51,71)=exp(y(71));
  g1(51,79)=(-(exp(y(62))*exp(y(79))))/((1+exp(y(79)))*(1+exp(y(79))));
  g1(51,84)=(-exp(y(84)));
  g1(51,101)=(-(exp(y(103))*exp(y(101))));
  g1(51,103)=(-(exp(y(103))*exp(y(101))));
  g1(52,47)=T1214;
  g1(52,50)=(-exp(y(50)));
  g1(52,5)=T730;
  g1(52,55)=T730;
  g1(52,62)=exp(y(62))*params(30);
  g1(52,9)=exp(y(9))*(1+exp(y(72)))/params(3);
  g1(52,63)=(-exp(y(63)));
  g1(52,64)=(-(exp(y(64))*exp(y(70))));
  g1(52,70)=(-(exp(y(64))*exp(y(70))));
  g1(52,72)=T746;
  g1(52,75)=(-exp(y(75)));
  g1(52,76)=(-exp(y(76)));
  g1(52,82)=1;
  g1(52,83)=(-(exp(y(83))*exp(y(84))));
  g1(52,84)=(-(exp(y(83))*exp(y(84))));
  g1(52,97)=1;
  g1(53,47)=(-((-((exp(y(64))*exp(y(70))-exp(y(62))*params(30))*T1219))/(T735*T735)));
  g1(53,62)=(-((-(exp(y(62))*params(30)))/T735));
  g1(53,64)=(-(exp(y(64))*exp(y(70))/T735));
  g1(53,70)=(-(exp(y(64))*exp(y(70))/T735));
  g1(53,71)=(-((-((exp(y(64))*exp(y(70))-exp(y(62))*params(30))*T735))/(T735*T735)));
  g1(53,98)=1;
  g1(54,47)=(-((T735*(-T1214)-T748*T1219)/(T735*T735)));
  g1(54,50)=(-(exp(y(50))/T735));
  g1(54,5)=(-((-T730)/T735));
  g1(54,55)=(-((-T730)/T735));
  g1(54,9)=(-(((-(exp(y(9))/params(3)))-T746)/T735));
  g1(54,63)=(-(exp(y(63))/T735));
  g1(54,71)=(-((-(T735*T748))/(T735*T735)));
  g1(54,72)=(-((-T746)/T735));
  g1(54,75)=(-(exp(y(75))/T735));
  g1(54,76)=(-(exp(y(76))/T735));
  g1(54,83)=(-(exp(y(83))*exp(y(84))/T735));
  g1(54,84)=(-(exp(y(83))*exp(y(84))/T735));
  g1(54,99)=1;
  g1(55,47)=(-(exp(y(47))*y(82)));
  g1(55,19)=(-(exp(y(19))*(1-params(51))/T118));
  g1(55,80)=exp(y(80))*(1-params(51));
  g1(55,81)=(-((-(T118*exp(y(24))))/(T118*T118)+(-(T118*exp(y(19))*(1-params(51))))/(T118*T118)-(-(T118*exp(y(27))*exp(y(93))))/(T118*T118)));
  g1(55,82)=(-exp(y(47)));
  g1(55,24)=(-(exp(y(24))/T118));
  g1(55,92)=exp(y(92));
  g1(55,93)=T314;
  g1(55,27)=T314;
  g1(55,95)=(-exp(y(95)));
  g1(56,42)=(-(exp(y(108))*T9*getPowerDeriv(T9,1-params(62),1)/(1-params(62))));
  g1(56,43)=(-(exp(y(116))*(-exp(y(43)))*getPowerDeriv(1-exp(y(43)),1-params(63),1)/(1-params(63))));
  g1(56,53)=(-(exp(y(114))*exp(y(53))*getPowerDeriv(exp(y(53)),1-params(65),1)/(1-params(65))));
  g1(56,80)=(-(exp(y(115))*exp(y(80))*getPowerDeriv(exp(y(80)),1-params(64),1)/(1-params(64))));
  g1(56,108)=(-T768);
  g1(56,114)=(-T777);
  g1(56,115)=(-T782);
  g1(56,116)=(-T772);
  g1(56,117)=exp(y(117));
  g1(56,138)=(-(params(2)*params(3)*exp(y(138))));
  g1(57,16)=(-(exp(x(it_, 1))*T796*exp(y(16))*getPowerDeriv(exp(y(16)),params(93),1)));
  g1(57,77)=exp(y(77));
  g1(57,139)=(-T800);
  g1(58,14)=(-(exp(x(it_, 2))*T808*exp(y(14))*getPowerDeriv(exp(y(14)),params(94),1)));
  g1(58,75)=exp(y(75));
  g1(58,140)=(-T812);
  g1(59,3)=(-(exp(x(it_, 3))*T820*exp(y(3))*getPowerDeriv(exp(y(3)),params(95),1)));
  g1(59,50)=exp(y(50));
  g1(59,141)=(-T824);
  g1(60,21)=(-(T832*exp(y(21))*getPowerDeriv(exp(y(21)),params(96),1)/exp(x(it_, 4))));
  g1(60,83)=exp(y(83));
  g1(60,142)=(-((-(T833*exp(x(it_, 4))))/(exp(x(it_, 4))*exp(x(it_, 4)))));
  g1(61,22)=(-(exp(x(it_, 5))*T844*exp(y(22))*getPowerDeriv(exp(y(22)),params(97),1)));
  g1(61,84)=exp(y(84));
  g1(61,143)=(-T848);
  g1(62,28)=(-(exp(x(it_, 6))*T856*exp(y(28))*getPowerDeriv(exp(y(28)),params(98),1)));
  g1(62,100)=exp(y(100));
  g1(62,144)=(-T860);
  g1(63,15)=(-(exp(x(it_, 7))*T868*exp(y(15))*getPowerDeriv(exp(y(15)),params(99),1)));
  g1(63,76)=exp(y(76));
  g1(63,145)=(-T872);
  g1(64,12)=(-(exp(x(it_, 9))*T880*exp(y(12))*getPowerDeriv(exp(y(12)),params(101),1)));
  g1(64,73)=exp(y(73));
  g1(64,147)=(-T884);
  g1(65,17)=(-(exp(x(it_, 10))*T892*exp(y(17))*getPowerDeriv(exp(y(17)),params(102),1)));
  g1(65,78)=exp(y(78));
  g1(65,148)=(-T896);
  g1(66,18)=(-(exp(x(it_, 11))*T904*exp(y(18))*getPowerDeriv(exp(y(18)),params(103),1)));
  g1(66,79)=exp(y(79));
  g1(66,149)=(-T908);
  g1(67,13)=(-(exp(x(it_, 12))*T916*exp(y(13))*getPowerDeriv(exp(y(13)),params(104),1)));
  g1(67,74)=exp(y(74));
  g1(67,150)=(-T920);
  g1(68,25)=(-(exp(x(it_, 13))*T928*exp(y(25))*getPowerDeriv(exp(y(25)),params(105),1)));
  g1(68,93)=exp(y(93));
  g1(68,151)=(-T932);
  g1(69,20)=(-(exp(x(it_, 14))*T939*getPowerDeriv(y(20),params(106),1)));
  g1(69,82)=1;
  g1(69,152)=(-T943);
  g1(70,10)=(-(exp(x(it_, 15))*T951*exp(y(10))*getPowerDeriv(exp(y(10)),params(107),1)));
  g1(70,69)=exp(y(69));
  g1(70,153)=(-T955);
  g1(71,29)=(-(T963*exp(y(29))*getPowerDeriv(exp(y(29)),params(108),1)/exp(x(it_, 16))));
  g1(71,107)=exp(y(107));
  g1(71,154)=(-((-(T964*exp(x(it_, 16))))/(exp(x(it_, 16))*exp(x(it_, 16)))));
  g1(72,30)=(-(T975*exp(y(30))*getPowerDeriv(exp(y(30)),params(109),1)/exp(x(it_, 17))));
  g1(72,108)=exp(y(108));
  g1(72,155)=(-((-(T976*exp(x(it_, 17))))/(exp(x(it_, 17))*exp(x(it_, 17)))));
  g1(73,31)=(-(exp(x(it_, 18))*T987*exp(y(31))*getPowerDeriv(exp(y(31)),params(110),1)));
  g1(73,109)=exp(y(109));
  g1(73,156)=(-T991);
  g1(74,34)=(-(exp(x(it_, 19))*T999*exp(y(34))*getPowerDeriv(exp(y(34)),params(111),1)));
  g1(74,112)=exp(y(112));
  g1(74,157)=(-T1003);
  g1(75,35)=(-(exp(x(it_, 20))*T1011*exp(y(35))*getPowerDeriv(exp(y(35)),params(112),1)));
  g1(75,113)=exp(y(113));
  g1(75,158)=(-T1015);
  g1(76,38)=(-(exp(x(it_, 21))*T1023*exp(y(38))*getPowerDeriv(exp(y(38)),params(113),1)));
  g1(76,116)=exp(y(116));
  g1(76,159)=(-T1027);
  g1(77,37)=(-(exp(x(it_, 22))*T1035*exp(y(37))*getPowerDeriv(exp(y(37)),params(114),1)));
  g1(77,115)=exp(y(115));
  g1(77,160)=(-T1039);
  g1(78,36)=(-(T1047*exp(y(36))*getPowerDeriv(exp(y(36)),params(115),1)/exp(x(it_, 23))));
  g1(78,114)=exp(y(114));
  g1(78,161)=(-((-(T1048*exp(x(it_, 23))))/(exp(x(it_, 23))*exp(x(it_, 23)))));
  g1(79,39)=(-(exp(x(it_, 24))*T1059*exp(y(39))*getPowerDeriv(exp(y(39)),params(116),1)));
  g1(79,120)=exp(y(120));
  g1(79,162)=(-T1063);
  g1(80,40)=(-(exp(x(it_, 25))*T1071*exp(y(40))*getPowerDeriv(exp(y(40)),params(117),1)));
  g1(80,121)=exp(y(121));
  g1(80,163)=(-T1075);
  g1(81,41)=(-(exp(x(it_, 26))*T1083*exp(y(41))*getPowerDeriv(exp(y(41)),params(118),1)));
  g1(81,122)=exp(y(122));
  g1(81,164)=(-T1087);
  g1(82,2)=params(45)*exp(y(9))*exp(y(2))/exp(y(11));
  g1(82,60)=exp(y(60));
  g1(82,9)=params(45)*exp(y(9))*exp(y(2))/exp(y(11));
  g1(82,11)=params(45)*(-(exp(y(9))*exp(y(2))*exp(y(11))))/(exp(y(11))*exp(y(11)));

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],82,26896);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],82,4410944);
end
end
end
end