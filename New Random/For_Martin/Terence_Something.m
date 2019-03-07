thick = 3;


x = 0:1:400;
line1 = 27504243 + 506546*x;

plot(x, line1, "LineWidth", thick)
hold on

cessna206 = 665000*x;
plot(x, cessna206, "LineWidth", thick)

cirrussr22TGTS = 696750*x;
plot(x, cirrussr22TGTS, "LineWidth", thick)


tecnum_p2006 = 597000*x;
plot(x, tecnum_p2006, "LineWidth", thick)

piper_PA46 = 917000*x;
plot(x, piper_PA46, "LineWidth", thick)

legend("Fixed + Variable Cost", "Cessna 206", "Cirrus SR22", "Tecnam P 2006", "Piper PA46")


text1 ="67";
text2 = "146";
text3 = "175";
text4 = "304";

x1 = 200;
y1 = 2000000;

x1 = 300;
y1 = 2000000;

x1 = 350;
y1 = 2000000;

x1 = 380;
y1 = 2000000;

% size = 20
% text(x1, y1, text1, "FontSize", size)
% text(x1, y1, text2, "FontSize", size)
% text(x1, y1, text3, "FontSize", size)
% text(x1, y1, text4, "FontSize", size)

xlabel("Numbers of Units Produced")
ylabel("Total Production Cost and Revenue")
set(gca, "FontSize", 30)