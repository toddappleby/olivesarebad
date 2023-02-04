%% loader
cd('/Users/toddappleby/Documents/Data/Clarinet Exports/SavedData')
load('MotionandNoise.mat')
%% grapher

figure(1);clf;hold on

%length
expCount = [];
for a = 1:length(MNONParasol.Exp)
    if ~isempty(MNONParasol.Exp(a).Cell)
        expCount = [expCount a];
    end
end
offset = .3;
for b = 1:length(expCount)
    
    subplot(2,round(length(expCount)/2),b)
    
    
      for bb = 1:length(MNONParasol.Exp(expCount(b)).Cell)
          plot(MNONParasol.Exp(expCount(b)).Cell(bb).LF.Seq(1:500)+(offset*(bb-1)),'Color','r')
          hold on
          
          plot(MNONParasol.Exp(expCount(b)).Cell(bb).LF.Rand(1:500)+(offset*(bb-1)),'Color','b')
          plot(MNONParasol.Exp(expCount(b)).Cell(bb).LF.Static(1:500)+(offset*(bb-1)),'Color','k')
          dateTitle = datestr(MNONParasol.Exp(expCount(b)).Cell(bb).Date);
          title(dateTitle(1:11))
      end
   set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.25, 9.125], 'PaperUnits', 'Inches', 'PaperSize', [5, 10],'Color','w')
      legend('Seq','Rand','Static','Location','northwest')
      set(gca,'xdir','reverse')
      sgtitle('ON Parasol Linear Filters')
      xlabel('time to spike (ms)')
      ylabel('weight')
end

figure(2);clf

offset = 100;


for b = 1:length(expCount)
    
    subplot(2,round(length(expCount)/2),b)
    
    
      for bb = 1:length(MNONParasol.Exp(expCount(b)).Cell)
          plot(MNONParasol.Exp(expCount(b)).Cell(bb).NL.X,(MNONParasol.Exp(expCount(b)).Cell(bb).NL.Seq*10000)+(offset*(bb-1)),'Color','r','LineWidth',2)
          hold on
          
          plot(MNONParasol.Exp(expCount(b)).Cell(bb).NL.X,(MNONParasol.Exp(expCount(b)).Cell(bb).NL.Rand)*10000+(offset*(bb-1)),'Color','b','LineWidth',2)
          plot(MNONParasol.Exp(expCount(b)).Cell(bb).NL.X,(MNONParasol.Exp(expCount(b)).Cell(bb).NL.Static)*10000+(offset*(bb-1)),'Color','k','LineWidth',2)
          dateTitle = datestr(MNONParasol.Exp(expCount(b)).Cell(bb).Date);
          title(dateTitle(1:11))
      end
   set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.25, 9.125], 'PaperUnits', 'Inches', 'PaperSize', [5, 10],'Color','w')
      legend('Seq','Rand','Static','Location','northwest')
%       set(gca,'xdir','reverse')
      sgtitle('ON Parasol Nonlinearities')
      xlabel('filtered resp')
      ylabel('spike rate (sp/s)')
end


figure(3);clf;hold on
expCount = [];
for a = 1:length(MNOFFParasol.Exp)
    if ~isempty(MNOFFParasol.Exp(a).Cell)
        expCount = [expCount a];
    end
end
offset = .3;
for b = 1:length(expCount)
    
    subplot(2,round(length(expCount)/2),b)
    
    
      for bb = 1:length(MNOFFParasol.Exp(expCount(b)).Cell)
          plot(MNOFFParasol.Exp(expCount(b)).Cell(bb).LF.Seq(1:500)+(offset*(bb-1)),'Color','r')
          hold on
          
          plot(MNOFFParasol.Exp(expCount(b)).Cell(bb).LF.Rand(1:500)+(offset*(bb-1)),'Color','b')
          plot(MNOFFParasol.Exp(expCount(b)).Cell(bb).LF.Static(1:500)+(offset*(bb-1)),'Color','k')
          dateTitle = datestr(MNOFFParasol.Exp(expCount(b)).Cell(bb).Date);
          title(dateTitle(1:11))
      end
   set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.25, 9.125], 'PaperUnits', 'Inches', 'PaperSize', [5, 10],'Color','w')
      legend('Seq','Rand','Static','Location','northwest')
      set(gca,'xdir','reverse')
      sgtitle('OFF Parasol Linear Filters')
      xlabel('time to spike (ms)')
      ylabel('weight')
end

figure(4);clf

offset = 100;


for b = 1:length(expCount)
    
    subplot(2,round(length(expCount)/2),b)
    
    
      for bb = 1:length(MNOFFParasol.Exp(expCount(b)).Cell)
          plot(MNOFFParasol.Exp(expCount(b)).Cell(bb).NL.X,(MNOFFParasol.Exp(expCount(b)).Cell(bb).NL.Seq*10000)+(offset*(bb-1)),'Color','r','LineWidth',2)
          hold on
          
          plot(MNOFFParasol.Exp(expCount(b)).Cell(bb).NL.X,(MNOFFParasol.Exp(expCount(b)).Cell(bb).NL.Rand)*10000+(offset*(bb-1)),'Color','b','LineWidth',2)
          plot(MNOFFParasol.Exp(expCount(b)).Cell(bb).NL.X,(MNOFFParasol.Exp(expCount(b)).Cell(bb).NL.Static)*10000+(offset*(bb-1)),'Color','k','LineWidth',2)
          dateTitle = datestr(MNOFFParasol.Exp(expCount(b)).Cell(bb).Date);
          title(dateTitle(1:11))
      end
   set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.25, 9.125], 'PaperUnits', 'Inches', 'PaperSize', [5, 10],'Color','w')
      legend('Seq','Rand','Static','Location','northwest')
%       set(gca,'xdir','reverse')
      sgtitle('OFF Parasol Nonlinearities')
      xlabel('filtered resp')
      ylabel('spike rate (sp/s)')
end



figure(5);clf
expCount = [];
for a = 1:length(MNONSmooth.Exp)
    if ~isempty(MNONSmooth.Exp(a).Cell)
        expCount = [expCount a];
    end
end
offset = .3;
for b = 1:length(expCount)
    
    subplot(2,round(length(expCount)/2),b)
    
    
      for bb = 1:length(MNONSmooth.Exp(expCount(b)).Cell)
          plot(MNONSmooth.Exp(expCount(b)).Cell(bb).LF.Seq(1:500)+(offset*(bb-1)),'Color','r')
          hold on
          
          plot(MNONSmooth.Exp(expCount(b)).Cell(bb).LF.Rand(1:500)+(offset*(bb-1)),'Color','b')
          plot(MNONSmooth.Exp(expCount(b)).Cell(bb).LF.Static(1:500)+(offset*(bb-1)),'Color','k')
          dateTitle = datestr(MNONSmooth.Exp(expCount(b)).Cell(bb).Date);
          title(dateTitle(1:11))
      end
   set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.25, 9.125], 'PaperUnits', 'Inches', 'PaperSize', [5, 10],'Color','w')
      legend('Seq','Rand','Static','Location','northwest')
      set(gca,'xdir','reverse')
      sgtitle('ON Smooth Linear Filters')
      xlabel('time to spike (ms)')
      ylabel('weight')
end

figure(6);clf

offset = 100;


for b = 1:length(expCount)
    
    subplot(2,round(length(expCount)/2),b)
    
    
      for bb = 1:length(MNONSmooth.Exp(expCount(b)).Cell)
          plot(MNONSmooth.Exp(expCount(b)).Cell(bb).NL.X,(MNONSmooth.Exp(expCount(b)).Cell(bb).NL.Seq*10000)+(offset*(bb-1)),'Color','r','LineWidth',2)
          hold on
          
          plot(MNONSmooth.Exp(expCount(b)).Cell(bb).NL.X,(MNONSmooth.Exp(expCount(b)).Cell(bb).NL.Rand)*10000+(offset*(bb-1)),'Color','b','LineWidth',2)
          plot(MNONSmooth.Exp(expCount(b)).Cell(bb).NL.X,(MNONSmooth.Exp(expCount(b)).Cell(bb).NL.Static)*10000+(offset*(bb-1)),'Color','k','LineWidth',2)
          dateTitle = datestr(MNONSmooth.Exp(expCount(b)).Cell(bb).Date);
          title(dateTitle(1:11))
      end
   set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.25, 9.125], 'PaperUnits', 'Inches', 'PaperSize', [5, 10],'Color','w')
      legend('Seq','Rand','Static','Location','northwest')
%       set(gca,'xdir','reverse')
      sgtitle('ON Smooth Nonlinearities')
      xlabel('filtered resp')
      ylabel('spike rate (sp/s)')
end



figure(7);clf
expCount = [];
for a = 1:length(MNOFFSmooth.Exp)
    if ~isempty(MNOFFSmooth.Exp(a).Cell)
        expCount = [expCount a];
    end
end
offset = .3;
for b = 1:length(expCount)
    
    subplot(2,round(length(expCount)/2),b)
    
    
      for bb = 1:length(MNOFFSmooth.Exp(expCount(b)).Cell)
          
          plot(MNOFFSmooth.Exp(expCount(b)).Cell(bb).LF.Seq(1:500)+(offset*(bb-1)),'Color','r')
          hold on
          
          plot(MNOFFSmooth.Exp(expCount(b)).Cell(bb).LF.Rand(1:500)+(offset*(bb-1)),'Color','b')
          plot(MNOFFSmooth.Exp(expCount(b)).Cell(bb).LF.Static(1:500)+(offset*(bb-1)),'Color','k')
          dateTitle = datestr(MNOFFSmooth.Exp(expCount(b)).Cell(bb).Date);
          title(dateTitle(1:11))
          
      end
   set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.25, 9.125], 'PaperUnits', 'Inches', 'PaperSize', [5, 10],'Color','w')
      legend('Seq','Rand','Static','Location','northwest')
      set(gca,'xdir','reverse')
      sgtitle('OFF Smooth Linear Filters')
      xlabel('time to spike (ms)')
      ylabel('weight')
end

figure(8);clf

offset = 100;


for b = 1:length(expCount)
    
    subplot(2,round(length(expCount)/2),b)
    
    
      for bb = 1:length(MNOFFSmooth.Exp(expCount(b)).Cell)
          plot(MNOFFSmooth.Exp(expCount(b)).Cell(bb).NL.X,(MNOFFSmooth.Exp(expCount(b)).Cell(bb).NL.Seq*10000)+(offset*(bb-1)),'Color','r','LineWidth',2)
          hold on
          
          plot(MNOFFSmooth.Exp(expCount(b)).Cell(bb).NL.X,(MNOFFSmooth.Exp(expCount(b)).Cell(bb).NL.Rand)*10000+(offset*(bb-1)),'Color','b','LineWidth',2)
          plot(MNOFFSmooth.Exp(expCount(b)).Cell(bb).NL.X,(MNOFFSmooth.Exp(expCount(b)).Cell(bb).NL.Static)*10000+(offset*(bb-1)),'Color','k','LineWidth',2)
          dateTitle = datestr(MNOFFSmooth.Exp(expCount(b)).Cell(bb).Date);
          title(dateTitle(1:11))
      end
   set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7.25, 9.125], 'PaperUnits', 'Inches', 'PaperSize', [5, 10],'Color','w')
      legend('Seq','Rand','Static','Location','northwest')
%       set(gca,'xdir','reverse')
      sgtitle('OFF Smooth Nonlinearities')
      xlabel('filtered resp')
      ylabel('spike rate (sp/s)')
end
%% MCS automation

cd('/Users/toddappleby/Documents/Data/Clarinet Exports/SavedData/CMS')
saveDir = dir;

    firstSlot=3;
    if strcmp(saveDir(3).name(1:3),'.DS')
        firstSlot = 4;
    end

    
    OFFPHolder = [];
    ONPHolder = [];
    ONSHolder = [];
    OFFSHolder = [] ;
    
    
for p = firstSlot:length(saveDir)

load(saveDir(p).name)
p
if ~isempty(CMSONParasol)
    for w = 1:size(CMSONParasol,1)
        ONPHolder(size(ONPHolder,1)+1,:) = CMSONParasol(w,:);
    end
end
if ~isempty(CMSOFFParasol)
    for x = 1:size(CMSOFFParasol,1)
        OFFPHolder(size(OFFPHolder,1)+1,:) = CMSOFFParasol(x,:);
    end
end
if ~isempty(CMSONSmooth)
    for y = 1:size(CMSONSmooth,1)
        
        ONSHolder(size(ONSHolder,1)+1,:) = CMSONSmooth(y,:);
    end
end
if ~isempty(CMSOFFSmooth)
    for z = 1:size(CMSOFFSmooth,1)
        OFFSHolder(size(OFFSHolder,1)+1,:) = CMSOFFSmooth(z,:);
    end
end

end

ONPHolder(isnan(ONPHolder))=0;
OFFPHolder(isnan(OFFPHolder))=0;
ONSHolder(isnan(ONSHolder))=0;
OFFSHolder(isnan(OFFSHolder))=0;

figure(10); clf;


meanArray = mean(OFFSHolder);
semArray = sem(OFFSHolder);

sgtitle('OFF Smooth Local/Global Motion')

subplot(2,3,1)
bar(1:2,[meanArray(1) meanArray(2)])
hold on
er = errorbar(1:2,[meanArray(1) meanArray(2)],[semArray(1) semArray(2)]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off
title('positive contrast center, separated surround')
ylabel('spike rate (sp/s)')
xlabel('1:sequential,2:random')

subplot(2,3,2)
bar(1:2,[meanArray(3) meanArray(4)])
hold on
er = errorbar(1:2,[meanArray(3) meanArray(4)],[semArray(3) semArray(4)]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off
title('negative contrast center, separated surround')
ylabel('spike rate (sp/s)')
xlabel('1:sequential,2:random')

subplot(2,3,3)
bar(1:3,[meanArray(5) meanArray(6) meanArray(7)])
hold on
er = errorbar(1:3,[meanArray(5) meanArray(6) meanArray(7)],[semArray(5) semArray(6) semArray(7)]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off
title('positive contrast separated center, sequential surround')
xlabel('1:seq,2:random,3:180')

subplot(2,3,4)
bar(1:3,[meanArray(8) meanArray(9) meanArray(10)])
hold on
er = errorbar(1:3,[meanArray(8) meanArray(9) meanArray(10)],[semArray(8) semArray(9) semArray(10)]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off
title('positive contrast separated center, random surround')
xlabel('1:seq,2:random,3:180')

subplot(2,3,5)
bar(1:3,[meanArray(11) meanArray(12) meanArray(13)])
hold on
er = errorbar(1:3,[meanArray(11) meanArray(12) meanArray(13)],[semArray(11) semArray(12) semArray(13)]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off
title('negtive contrast separated center, sequential surround')
xlabel('1:seq,2:random,3:180')

subplot(2,3,6)
bar(1:3,[meanArray(14) meanArray(15) meanArray(16)])
hold on
er = errorbar(1:3,[meanArray(14) meanArray(15) meanArray(16)],[semArray(14) semArray(15) semArray(16)]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off
title('negative contrast separated center, random surround')
xlabel('1:seq,2:random,3:180')

figure(11); clf;

meanArray = mean(ONSHolder);
semArray = sem(ONSHolder);

sgtitle('ON Smooth Local/Global Motion')

subplot(2,3,1)
bar(1:2,[meanArray(1) meanArray(2)])
hold on
er = errorbar(1:2,[meanArray(1) meanArray(2)],[semArray(1) semArray(2)]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off
title('positive contrast center, separated surround')
ylabel('spike rate (sp/s)')
xlabel('1:sequential,2:random')

subplot(2,3,2)
bar(1:2,[meanArray(3) meanArray(4)])
hold on
er = errorbar(1:2,[meanArray(3) meanArray(4)],[semArray(3) semArray(4)]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off
title('negative contrast center, separated surround')
xlabel('1:sequential,2:random')

subplot(2,3,3)
bar(1:3,[meanArray(5) meanArray(6) meanArray(7)])
hold on
er = errorbar(1:3,[meanArray(5) meanArray(6) meanArray(7)],[semArray(5) semArray(6) semArray(7)]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off
title('positive contrast separated center, sequential surround')
xlabel('1:seq,2:random,3:180')

subplot(2,3,4)
bar(1:3,[meanArray(8) meanArray(9) meanArray(10)])
hold on
er = errorbar(1:3,[meanArray(8) meanArray(9) meanArray(10)],[semArray(8) semArray(9) semArray(10)]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off
title('positive contrast separated center, random surround')
xlabel('1:seq,2:random,3:180')

subplot(2,3,5)
bar(1:3,[meanArray(11) meanArray(12) meanArray(13)])
hold on
er = errorbar(1:3,[meanArray(11) meanArray(12) meanArray(13)],[semArray(11) semArray(12) semArray(13)]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off
title('negtive contrast separated center, sequential surround')
xlabel('1:seq,2:random,3:180')

subplot(2,3,6)
bar(1:3,[meanArray(14) meanArray(15) meanArray(16)])
hold on
er = errorbar(1:3,[meanArray(14) meanArray(15) meanArray(16)],[semArray(14) semArray(15) semArray(16)]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off
title('negative contrast separated center, random surround')
xlabel('1:seq,2:random,3:180')

figure(12); clf;

meanArray = mean(OFFPHolder);
semArray = sem(OFFPHolder);

sgtitle('OFF Parasol Local/Global Motion')

subplot(2,3,1)
bar(1:2,[meanArray(1) meanArray(2)])
hold on
er = errorbar(1:2,[meanArray(1) meanArray(2)],[semArray(1) semArray(2)]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off
title('positive contrast center, separated surround')
ylabel('spike rate (sp/s)')
xlabel('1:sequential,2:random')

subplot(2,3,2)
bar(1:2,[meanArray(3) meanArray(4)])
hold on
er = errorbar(1:2,[meanArray(3) meanArray(4)],[semArray(3) semArray(4)]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off
title('negative contrast center, separated surround')
xlabel('1:sequential,2:random')

subplot(2,3,3)
bar(1:3,[meanArray(5) meanArray(6) meanArray(7)])
hold on
er = errorbar(1:3,[meanArray(5) meanArray(6) meanArray(7)],[semArray(5) semArray(6) semArray(7)]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off
title('positive contrast separated center, sequential surround')
xlabel('1:seq,2:random,3:180')

subplot(2,3,4)
bar(1:3,[meanArray(8) meanArray(9) meanArray(10)])
hold on
er = errorbar(1:3,[meanArray(8) meanArray(9) meanArray(10)],[semArray(8) semArray(9) semArray(10)]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off
title('positive contrast separated center, random surround')
xlabel('1:seq,2:random,3:180')

subplot(2,3,5)
bar(1:3,[meanArray(11) meanArray(12) meanArray(13)])
hold on
er = errorbar(1:3,[meanArray(11) meanArray(12) meanArray(13)],[semArray(11) semArray(12) semArray(13)]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off
title('negtive contrast separated center, sequential surround')
xlabel('1:seq,2:random,3:180')

subplot(2,3,6)
bar(1:3,[meanArray(14) meanArray(15) meanArray(16)])
hold on
er = errorbar(1:3,[meanArray(14) meanArray(15) meanArray(16)],[semArray(14) semArray(15) semArray(16)]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off
title('negative contrast separated center, random surround')
xlabel('1:seq,2:random,3:180')

figure(13); clf;

meanArray = mean(ONPHolder);
semArray = sem(ONPHolder);

sgtitle('ON Parasol Local/Global Motion')

subplot(2,3,1)
bar(1:2,[meanArray(1) meanArray(2)])
hold on
er = errorbar(1:2,[meanArray(1) meanArray(2)],[semArray(1) semArray(2)]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off
title('positive contrast center, separated surround')
ylabel('spike rate (sp/s)')
xlabel('1:sequential,2:random')

subplot(2,3,2)
bar(1:2,[meanArray(3) meanArray(4)])
hold on
er = errorbar(1:2,[meanArray(3) meanArray(4)],[semArray(3) semArray(4)]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off
title('negative contrast center, separated surround')
xlabel('1:sequential,2:random')

subplot(2,3,3)
bar(1:3,[meanArray(5) meanArray(6) meanArray(7)])
hold on
er = errorbar(1:3,[meanArray(5) meanArray(6) meanArray(7)],[semArray(5) semArray(6) semArray(7)]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off
title('positive contrast separated center, sequential surround')
xlabel('1:seq,2:random,3:180')

subplot(2,3,4)
bar(1:3,[meanArray(8) meanArray(9) meanArray(10)])
hold on
er = errorbar(1:3,[meanArray(8) meanArray(9) meanArray(10)],[semArray(8) semArray(9) semArray(10)]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off
title('positive contrast separated center, random surround')
xlabel('1:seq,2:random,3:180')

subplot(2,3,5)
bar(1:3,[meanArray(11) meanArray(12) meanArray(13)])
hold on
er = errorbar(1:3,[meanArray(11) meanArray(12) meanArray(13)],[semArray(11) semArray(12) semArray(13)]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off
title('negtive contrast separated center, sequential surround')
xlabel('1:seq,2:random,3:180')

subplot(2,3,6)
bar(1:3,[meanArray(14) meanArray(15) meanArray(16)])
hold on
er = errorbar(1:3,[meanArray(14) meanArray(15) meanArray(16)],[semArray(14) semArray(15) semArray(16)]);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off
title('negative contrast separated center, random surround')
xlabel('1:seq,2:random,3:180')
