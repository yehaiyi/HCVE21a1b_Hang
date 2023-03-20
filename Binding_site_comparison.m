%% 获取h77和consensus 1b
load("data_new\h77.mat")
load("data_new/Cseq_mut_order_1a&1b.mat")


consensus_1b = Cseq_mut_order_1b ;

%%
% 获取 consensus 3a
load("data_new\NumofPatient_3aE2.mat")
inputfile = 'data_new/3a_E2_ori.fasta';
[Header_fasta, Sequence_fasta] = fastaread(inputfile);
msa_aa = cell2mat(Sequence_fasta');

% preprocess 
no_patient_idx = find(~patient);
if length(no_patient_idx) >0
    msa_aa(no_patient_idx) = []
end

load("data_new\outliers.mat")   % btw, 这里的outliers和NumofPatient里面带的不一样，所以是咋来的来着
msa_aa(outliers,:) =[];

num_res= size(msa_aa,2);
unique_aa =[];
ind_conserve = [];
ind_non_conserve=[];
consensus_3a = [];
for i = 1:num_res
    col = msa_aa(:,i);
    unique_aa = [unique_aa; string(unique(col)') ];
    if length(unique(col)) == 1
        ind_conserve = [ind_conserve i]; 
    else 
        ind_non_conserve=[ind_non_conserve i];
    end
    consensus_3a = [consensus_3a '1'];
end

load("data_new\3aCseq_mut_order.mat")
consensus_3a(ind_non_conserve)=Cseq_mut_order;
consensus_3a(ind_conserve) = unique_aa(ind_conserve)


%% Binding sites for 35 Abs

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

HmAbs_all=[HmAbs_all  HmAbs_Bailey2017 HmAbs_Gopal2017 HmAbs_Keck2019 HmAb_Pierce2016];

%% 
load("data_new\seqs_h77_cons1b_cons3a.mat")
load("data_new\McLachlan_similarity.mat")
Diff_res_all = [];
Same_res_all =[];
sites_all =[];

for cell = 1:35
    
    sites = HmAbs_all{1,cell};
    sites_all = [sites_all sites];
    sites_3a=[];
    % index转换
    for site = sites
       
        % convert H77 number to S52 number
        if site>479 && site<579
            site=site+1;
            sites_3a = [sites_3a site];
        else if site>578
            site=site+6;
            sites_3a = [sites_3a site];
        else
            sites_3a= [sites_3a site];
        end
    end
    end 
  

    bs_for_3_strains = [h77(sites -383) ; consensus_1b(sites - 383) ; consensus_3a(sites_3a -383)];

    % print aa of these binding sites
    for i=1:3
        show_aa =[];
        for j = bs_for_3_strains(i,:)
            show_aa = [show_aa ',' j];  
        end
        show_aa;
    end

    % count inconsistent residues
    same_count =0;
    diff_res = [];
    for i = 1: length(sites)
        if length(unique(bs_for_3_strains(:,i))) ==1
            same_count =  same_count+ 1;
            Same_res_all = [Same_res_all sites(i)];
        else
            diff_res = [diff_res sites(i)];
        end
    end 
    Diff_res_all = [ Diff_res_all diff_res];  %这里是以h77为ref的number system

    Num_2 = length(find((bs_for_3_strains(1,:) ~= bs_for_3_strains(2,:))));
    Num_3 =  length(sites)- same_count;

    % caculate chemical similarity for binding sites of each antibody
    s1a1b =[];
    s1a3a = [];
    s1b3a = [];
    for j = 1: size(bs_for_3_strains,2)
        pair1a1b = [ bs_for_3_strains(1,j) bs_for_3_strains(2,j) ]; 
        pair1a3a = [ bs_for_3_strains(1,j) bs_for_3_strains(3,j) ]; 
        pair1b3a = [bs_for_3_strains(2,j) bs_for_3_strains(3,j)];
        if isfield(scores,pair1a1b)
            s = getfield(scores,pair1a1b);
            %s1a1b = [s1a1b ',' sprintf("%d",s)];
            s1a1b = [ s1a1b s];
        else
            s=0;
            s1a1b = [s1a1b s];
        end
        if isfield(scores,pair1a3a)
            s = getfield(scores,pair1a3a);
            s1a3a = [s1a3a s ];
        else 
            s =0;
            s1a3a = [s1a3a s];
        end
        if isfield(scores,pair1b3a)
            s = getfield(scores,pair1b3a);
            s1b3a = [s1b3a s ];
       else 
            s =0;
            s1b3a = [s1b3a s];
        end

    end

    sprintf("cell=%d" , cell)
    average_CSS = [mean(s1a1b)  mean(s1a3a)  mean(s1b3a)   mean([mean(s1a1b)  mean(s1a3a)  mean(s1b3a) ])]
    
    %s1a1b
    %s1a3a
    %sprintf(" %d/%d, %d/%d ", Num_2,length(sites), Num_3, length(sites))
    %sprintf("length:=%d same=%d," , length(sites), same_count)   %
    %inconsistent

end

Diff_res_all = unique(Diff_res_all);
Diff_res_all_3a=[];
for site = Diff_res_all;
       
        % convert H77 number to S52 number
        if site>479 && site<579
            site=site+1;
            Diff_res_all_3a = [Diff_res_all_3a site];
        else if site>578
            site=site+6;
            Diff_res_all_3a = [Diff_res_all_3a site];
        else
            Diff_res_all_3a= [Diff_res_all_3a site];
        end
    end
end 

% show inconsisitent aa
Diff_aa = [h77(Diff_res_all -383)  ; consensus_1b(Diff_res_all - 383 )  ;  consensus_3a(Diff_res_all_3a - 383)] ;

for i=1:3
        show_aa =[];
        for j = Diff_aa(i,:)
            show_aa = [show_aa ',' j];
        end
        show_aa;
end

sites_all = unique(sites_all);
Same_res_all = unique(Same_res_all);

%%  chemical similarity
% 对涉及到的所有 binding sites (127个residues) / inconsistent binding sites
%sites_all  = unique(Diff_res_all);
s1a1b =[];
s1a3a = [];
s1b3a = [];

sites_all_3a = [];
for site = sites_all;
       % convert H77 number to S52 number
        if site>479 && site<579
            site=site+1;
            sites_all_3a = [sites_all_3a site];
        else if site>578
            site=site+6;
            sites_all_3a = [sites_all_3a site];
        else
            sites_all_3a= [sites_all_3a site];
        end
    end
end 

for k = 1 : length(sites_all)
     j = sites_all(k) - 383 ;
     j_3a = sites_all_3a(k)- 383;
     pair1a1b = [ h77(j) consensus_1b(j) ]; 
     pair1a3a = [ h77(j) consensus_3a(j_3a) ]; 
     pair1b3a = [consensus_1b(j) consensus_3a(j_3a)];
     if isfield(scores,pair1a1b)
         s = getfield(scores,pair1a1b);
         %s1a1b = [s1a1b ',' sprintf("%d",s)];
         s1a1b = [ s1a1b s];
     else
         s=0;
         s1a1b = [s1a1b s];
     end
     if isfield(scores,pair1a3a)
         s = getfield(scores,pair1a3a);
         s1a3a = [s1a3a s ];
     else 
         s =0;
         s1a3a = [s1a3a s];
     end
     if isfield(scores,pair1b3a)
         s = getfield(scores,pair1b3a);
         s1b3a = [s1b3a s ];
     else 
         s =0;
         s1b3a = [s1b3a s];
     end
end
s1b3a

%% 将binding residues进一步分为escape-resistent和not escape-resistent
load("data_new\escape-resistent_residues_both.mat")

re =[];
nre = [];

for site = sites_all
    if resistent_both(site-383)==1
        re = [re site];
    else 
        nre = [nre site];
    end
end

length(re)
length(nre)
%% escape-resistent binding sites 和 not escape-resistent binding sites的average CCS

sites_all =Same_res_all ;  % re or nre

s1a1b =[];
s1a3a = [];
s1b3a = [];

sites_all_3a = [];
for site = sites_all;
       % convert H77 number to S52 number
        if site>479 && site<579
            site=site+1;
            sites_all_3a = [sites_all_3a site];
        else if site>578
            site=site+6;
            sites_all_3a = [sites_all_3a site];
        else
            sites_all_3a= [sites_all_3a site];
        end
    end
end 

for k = 1 : length(sites_all)
     j = sites_all(k) - 383 ;
     j_3a = sites_all_3a(k)- 383;
     pair1a1b = [ h77(j) consensus_1b(j) ]; 
     pair1a3a = [ h77(j) consensus_3a(j_3a) ]; 
     pair1b3a = [consensus_1b(j) consensus_3a(j_3a)];
     if isfield(scores,pair1a1b)
         s = getfield(scores,pair1a1b);
         %s1a1b = [s1a1b ',' sprintf("%d",s)];
         s1a1b = [ s1a1b s];
     else
         s=0;
         s1a1b = [s1a1b s];
     end
     if isfield(scores,pair1a3a)
         s = getfield(scores,pair1a3a);
         s1a3a = [s1a3a s ];
     else 
         s =0;
         s1a3a = [s1a3a s];
     end
     if isfield(scores,pair1b3a)
         s = getfield(scores,pair1b3a);
         s1b3a = [s1b3a s ];
     else 
         s =0;
         s1b3a = [s1b3a s];
     end
end

length(s1a1b)
mean(s1a1b)
mean(s1a3a)
mean(s1b3a)

%% 
load("data_new\escape_time.mat")

time_1a = escape_time(1,:);
time_1b = escape_time(2,:);

%% 
escape_residues_1a=[384   386   388   390   391   393   394   395   396   397   398   399   400 ...
                    401   402   403   404   405   407   408   410   415   416   417   422   424 ...
                    431   433   434   435   438   442   444   446   453   456   461   466   475 ...
                    482   501   524   528   531   533   538   557   558   560   580   608 ...
                    610   636   713];  % table + 439 578
%% 
poly=[416, 422, 424, 431, 433, 438, 442, 446, 453, 456, 461, 475, 482, 520 , 524, 531, 533, 557, 558, 560]

%% exclude Domain A antibodies
domain_a = [626:632];

for cell = 1:35
    binding_sites = HmAbs_all{1,cell};
    inter = intersect(domain_a, binding_sites);
    if isempty(inter)
        sprintf("cell=%d, NAbs", cell)
    elseif length(inter) == length(binding_sites)
        sprintf("cell=%d, non NAbs", cell)
    else
        sprintf("cell=%d, some in Domain A", cell)
    end

end
