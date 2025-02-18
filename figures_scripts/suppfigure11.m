HmAbs_all=[];

HmAb_Pierce2016 = [];

%CBH-4D
HmAb_Pierce2016{1} = [494 497 502 504 506:509 511 537 539 540 542:545 547 549:552 554 556 559 561 562 564 565 584 585 592 594 598 600 602 603 607:611 614 617:619 621 623 624 626 627 629:633 638 640 642:644];
%CBH-4G
HmAb_Pierce2016{2} = [494 497 503 505:509 511 517 537 539 540 542:545 547 549:552 554 556 559 561 564 565 584 592:594 597 598 600 602 603 607:611 614 617:619 621 623 624 627 629 631:633 638 640 642:644];
%CBH-4B
HmAb_Pierce2016{3} = [494 497 503 505:509 537 539 540 550:552 554 559 561 564 565 600 602 603 607:611 614 617:619 621 624 627 629 631:633 638 640 642:644];
%CBH-20
% HmAb_Pierce2016{4} = [494 497 503 505:509 537 539 540 542:545 547 549:552 554 559 561 564 565 569 581 583:585 592 594 597 598 600 602 603 607 608 611 614 617:619 621 623 624 627 631 632 638 640 642:644];
%included the residues with close to threshold (RB = 21,22) as well
HmAb_Pierce2016{4} = [494 497 503 505:509 537 539 540 542:545 547 549:552 554 559 561 564 565 569 581 583:585 592 594 597 598 600 602 603 607 608 610 611 614 617:619 621 623 624 627 631 632 638 640 642:644]; 
%CBH-21
HmAb_Pierce2016{5} = [452 494 497 503 505:509 537 539 540 542:545 547 549:552 554 559 561 564 565 569 581 583:585 592:594 597 598 600 602 603 607:611 614 617:619 621 623 624 627 631:633 638 640 642:644];
%CBH-22
HmAb_Pierce2016{6} = [452 494 497 503 505:509 537 539 540 542:545 547 549:552 554 559 561 564 565 569 581 583:585 592:594 597 598 600 602 603 607:611 614 617:619 621 623 624 627 631:633 638 640 642:644];

%HC-1




HmAb_Pierce2016{7} = [429 494 503:506 508 509 529 530 535 537 539 552 554 559 564 607 611 614 617 644];
%HC-11
HmAb_Pierce2016{8} = [425 428 429 436:438 442 443 494 497 502:504 506:509 511 520 530 535 537 539 550:552 554 556 559 564 565 602 603 607 608 611 614 617:621 624 640 643 644];
%A27
HmAb_Pierce2016{9} = [424:429 437 438 494 497 499 502:504 506:509 511 520 529 530 535 537 539 540 550:552 554 556 559 564 565 602 603 607:611 614 616:619 621 623 624 638 640 642:644];

%CBH-23
HmAb_Pierce2016{10} = [494 508 509 537 539 549 552 554 564 611 614 644]; %same as CBH7
%CBH-7
HmAb_Pierce2016{11} = [494 506 508 509 537 539 549 552 554 564 611 614 617 621 644];

%HC84-20
HmAb_Pierce2016{12} = [429 441 494 497 502:509 511 537 539 552 554 559 564 607 608 611 613 614 616:619 621 640 643 644];
%HC84-24
HmAb_Pierce2016{13} = [429 442 443 494 497 502:509 511 537 539 552 554 559 564 607 608 611 614 617:619 621 643 644];
%HC84-26
HmAb_Pierce2016{14} = [441 442 494 497 502:509 511 537 539 552 554 559 564 603 607 608 611 614 617:619 621 640 643 644];

%HC33-1
HmAb_Pierce2016{15} = [413 418 420];
%HC33-4
HmAb_Pierce2016{16} = [408 413 420];
% %CD81_bs
% HmAb_Pierce2016{17} = [420 421 424 427 430 436:438 440:443 523 527 529 530 535 540 613 614 616:618];



%HCV1
HmAbs_Gopal2017{1} = [413 415 418 420 422 624];
%AR1A
HmAbs_Gopal2017{2} = [417 485 494 497 502 504 506:511 514 516 519 537:540 542 544 545 547 549:552 554 559 564 565 602 603 607:611 614 617:619 621 623 624 632 638 640 642 643 644]; %559 is RB=21%
%AR1B
HmAbs_Gopal2017{3} = [494 497 504 506:509 511 537 539 544 545 547:552 554 564 607:608 610 611 614 617:619 621 640 644]; %610 is RB=21%
%AR2A
HmAbs_Gopal2017{4} = [552 607:611 614 617:619 621 623:625 628 638 640 643 644];
%AR3A
HmAbs_Gopal2017{5} = [424 425 428 429 436:438 440:442 485 494 496 497 502:509 511 516 523 525 529 530 535 537 539 552 554 559 564 602 603 607:611 614 616:619 621 623 624 625 640 642:644];
%AR3B
HmAbs_Gopal2017{6} = [424 425 427:429 432 436 437 440:442 494 496 497 502:509 511 516 518 520 523 529 530 535 537 539 540 552 554 555 556 558 559 564 602 603 607:611 614 616:619 621 623 624 625 630 638 640:644];
%AR3C
HmAbs_Gopal2017{7} = [424 425 428:429 437 438 441:443 494 496 497 502:509 511 516 525 529 530 535 537 539 540 552 554 559 564 602 603 607:611 614 616:619 621 623 624 625 630 638 640 641 643 644];
%AR3D
HmAbs_Gopal2017{8} = [424 425 427:430 436:438 440:442 459 494 496 497 499 502:509 511 516 518 520 523 529 530 535 537 539 540 550:552 554 555 556 558 559 564 565 600 602 603 607:612 614 616:619 621 623 624 625 638 640:644];
 
 

%CBH-5
HmAbs_Keck2019{1} = [424:429 436:438 441:443 535 638];
%212.1.1
HmAbs_Keck2019{2} = [425:429 433 434 529 530 535];
%212.1.10
HmAbs_Keck2019{3} = [426 428 429 442];
%212.15
HmAbs_Keck2019{4} = [544 547 549 637 638 639];
%212.25
HmAbs_Keck2019{5} = [549 636 638 639];

 
% Bailey2017

%HEPC3
HmAbs_Bailey2017{1} = [425 427 428 437 499 520 530 535];
%HEPC43
HmAbs_Bailey2017{2} = [425 427 428 432 436 437 438 442 443 499 517 520 527 529 530 535 616];
%HEPC74
HmAbs_Bailey2017{3} = [425 428 436 437 530 535];
%HEPC46
HmAbs_Bailey2017{4} = [541:546 548 549 594 598 633];
%HEPC50
HmAbs_Bailey2017{5} = [543 544 545 549 594 597 598];
%HEPC98
HmAbs_Bailey2017{6} = [402 405 408];

HmAbs_all=[HmAbs_all HmAbs_Bailey2017  HmAbs_Gopal2017 HmAbs_Keck2019 HmAb_Pierce2016];

% HmAbs_all=[HmAbs_all HmAbs_Bailey2017  HmAbs_Gopal2017 HmAbs_Keck2019 HmAb_Pierce2016];
% 



label_xaxis_data =  {'HEPC3',...
    'HEPC43',' HEPC74',' HEPC46','HEPC50','HEPC98','HCV1','AR1A','AR1B','AR2A','AR3A','AR3B',...
    'AR3C','AR3D','CBH-5','212.1.1',' 212.10','212.15','212.25','CBH-4D','CBH-4G','CBH-4B','CBH-20','CBH-21','CBH-22',...
    'HC-1','HC-11',' A27','CBH-23','CBH-7',' HC84-20','HC84-24','HC84-26',...
    'HC33-1','HC33-4'};

label_xaxis_data =flip(label_xaxis_data);
HmAbs_all = flip(HmAbs_all);

label_yaxis_data = {'Subtype 1a','Subtype 1b'};

run startup.m
load('escape_time_1a_500.mat', 'mean_escape_time')
all_E=mean_escape_time;

map = zeros(length(label_yaxis_data),length(label_xaxis_data));
for i =1:length(label_yaxis_data)
for kk = 1:length(HmAbs_all)
    data = all_E(HmAbs_all{kk}-383);
%     map(i,kk) = length(find(data<= prctile(all_E,i*5)));
    map(1,kk) = length(find(data<=100));
end
end
load('mean_escape_time_1b.mat')
all_E=all_mean_escape_time_1b;
for kk = 1:length(HmAbs_all)
    data = all_E(HmAbs_all{kk}-383);
%     map(i,kk) = length(find(data<= prctile(all_E,i*5)));
    map(2,kk) = length(find(data<=80));
end




% map  = flip(map,1);
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)

FIG=figure;
% subplot(2,1,2)
[grad,im]=colorGradient(	hex2rgb('#FCFAF2'),hex2rgb('#4A225D'),128);
heatmap(label_xaxis_data,label_yaxis_data,map,'GridVisible','off','Colormap',grad,'ColorLimits',[0 max(max(map))],'ColorbarVisible','off','FontName','Arial','FontSize',8);
% text('Units', 'Normalized', 'Position', [-0.09,0.5],'String','Escape time threshold, \tau','Rotation',90)
% ,'FontName','Arial','FontSize',8,'FontWeight','Bold'

FIG.Units = 'centimeters';
FIG.Name = 'RBthresh';
set(gcf,'Position',[10 10 20 3]);
% ylabel('Escape time threshold, \tau')
set(gca,'Position',[.09 .4 .88 .5]);  %调整 XLABLE和YLABLE不会被切掉
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% print(['C:\Users\hzhangbr\Desktop\' FIG.Name],'-dpng','-r600');


FIG=figure;
FIG.Units = 'centimeters';
FIG.Name = 'xlabel';
xtickangle(45)
xlim([0.5 35.5])
set(gcf,'Position',[10 10 20 3]);
% ylabel('Escape time threshold, \tau')
set(gca,'Position',[.09 .4 .88 .5]);  %调整 XLABLE和YLABLE不会被切掉

xtick =set(gca,'XTick',1:35,'XTickLabel',...
    flip({'\bf HEPC3',...
    'HEPC43','\bf HEPC74','\bf HEPC46','HEPC50','HEPC98','HCV1','AR1A','AR1B','AR2A','AR3A','AR3B',...
    'AR3C','AR3D','CBH-5','212.1.1','\bf 212.10','212.15','212.25','CBH-4D','CBH-4G','CBH-4B','CBH-20','CBH-21','CBH-22',...
    'HC-1','HC-11',' A27','CBH-23','CBH-7','\bf HC84-20','HC84-24','HC84-26',...
    'HC33-1','HC33-4'}),'FontName','Arial');
% print(['C:\Users\hzhangbr\Desktop\' FIG.Name],'-dpng','-r600');
% print(['C:\Users\hzhangbr\Desktop\' FIG.Name],'-dpng','-r600');
% print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600');