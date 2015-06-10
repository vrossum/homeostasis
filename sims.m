%clear all
close all

% determinants of filter cascades

source('funs.m')
source('simfuns.m')

global glob_functionplot

iplot=1;
res=[];
simcase=31


switch simcase
  case{0} % dummy to re-read functions
	1;
 case{1} % effect of noise
	t_ar=[ 10, 50, 100] 
%	t_ar=[ 1, 2 ,5 ]
	nn=10
	ww=rand(nn);
	ww=ww+ww';
	maxlam=max(eig(ww))
	wwbase = ww./maxlam;
	resall=[];
	for si= [0.001 0.1]
	si
	res=[];
	tend=20000;

	 wscalear=[0.65:0.005:0.85];
	 for wscale=wscalear
		r1=sim_dyn_noisefull(wscale*wwbase,t_ar,tend,si);
		n2=length(r1)/2;
		osc_ampl=mean(std(r1(:,n2:end)'))
		res=[res [wscale osc_ampl]'];
	  end
	  figure(2)
	  plot(res(1,:),res(2,:),'x-');
	  hold on
	  title('osc_ampl vs w')
	  
	  % compare to without homeo
	  %note, above critical w =1, 1D system satuarates so no lvucs
	  compareQ=1;
	  eps=1e-2;
	  if (compareQ && si > eps)
		res2=[];
		for wscale=wscalear
		  r1=sim_dyn_noisefull(wscale*wwbase,t_ar(1),tend,si);
		  n2=round(length(r1)/2);
		  osc_ampl=mean(std(r1(:,n2:end)'))
		  res2=[res2 [wscale osc_ampl]'];
		end
		figure(2)
		plot(res2(1,:),res2(2,:),'r;non-homeo;');
     	  resall=[resall; res2'];
	  end 
	  resall=[resall; res']
	end
	save noise.dat resall
 
 case{2} % heterogeneity
	warning ("off", "Octave:broadcast");
	nn=10
	ww=rand(nn);
	ww=ww+ww';
	maxlam=max(eig(ww))
	wwbase = ww./maxlam;

	
	t_ar1d=[1 2 5];
	wc=get_wcrit(t_ar1d)
	
	ipl=1;
	for cv_tau=[0.01 0.1 0.2 0.3 0.4 0.5 0.6]
	res=[];
	for wscale=0.5:0.02:0.99
	  ww=wwbase*wscale;
	  kk=length(t_ar1d);
	  for i=1:kk 
	  % gamrnd(k,b): mu = k*b, var=k*b^2, cv=1/sqrt(k) 
		k = 1/cv_tau^2;
		t_ar2d(i,:)=gamrnd(k,t_ar1d(i)/k,1,nn);
	  end
	  %mmf=create_mmatfullhet(ww,t_ar2d)
	  tend=1000;
	  r1mat 	= fullsim_dyn_het(ww,t_ar2d,tend);
	  r1=r1mat(1,:);
	  n2		= round(length(r1)/2);
	  osc_ampl= std(r1(n2:end));
	  res=[res [wscale osc_ampl]'];
	end
	figure(2)
	pltstr=strcat(num2str(ipl),';',num2str(cv_tau),';'); % gives unique color line + legend.
	ipl++
	plot(res(1,:),res(2,:),pltstr);
	hold on
	end
	
 case{21} % heterogeneity based on eigenvectors
	warning ("off", "Octave:broadcast");
	imagQ=0; % criterion for oscilation
	nn=10
	ww=randn(nn);
	ww=ww+ww';
	maxlam=max(eig(ww))
	wwbase = ww./maxlam;
		
	t_ar1d=[10 50 100];
	%t_ar1d=[1 2 5];
	wc=get_wcrit(t_ar1d)
	
	fp=fopen("het_dat.dat","w");
	
	ipl=1;
	res=[];
	for cv_tau=0.01:0.01:1
	  k=1/cv_tau^2;
	  ntr=100;
	  wc_ar=[];
	  for itr=1:ntr
	  % could redraw ww here.
		kk=length(t_ar1d);
		for i=1:kk  % gamrnd(k,b): mu = k*b, var=k*b^2, cv=1/sqrt(k) 
		  t_ar2d(i,:)=gamrnd(k,t_ar1d(i)/k,1,nn);
		end
		wsmin=0; wsmax=1;
		myFunc1=@(wscale)max(real(eig(create_mmatfullhet(wwbase*wscale,t_ar2d))));
		
		if (imagQ) 
		  myFunc1=@(wscale)(-eps+max(abs(imag(eig(create_mmatfullhet(wwbase*wscale,t_ar2d))))));
		end
		if (sign(myFunc1(wsmin)*myFunc1(wsmax))>0)
		  fprintf("warning no bracket\n");
		  wc = (myFunc1(wsmax)<0);
		else
		  wc= fzero(myFunc1, [wsmin, wsmax]);
		end
		wc_ar = [wc_ar wc];	  
		fprintf(fp,"%g %g\n",cv_tau, wc);
	  end
	  plot(cv_tau,wc_ar,"*r")
	  drawnow();
	  hold on
	  res=[res [cv_tau mean(wc_ar) std(wc_ar) min(wc_ar)]'];
	end
	errorbar(res(1,:),res(2,:),res(3,:))
	xlabel("cv_tau");
	ylabel("mean wc");
	fclose(fp);
	rest=res';
	save het.dat rest
	
case{3} % non-linear
	t_ar=[ 10, 50, 10]
	
	  
	res=[];
	alpha=1;
	xthr=0.1;
	w=0;
	
	t3c=t_ar(1)*t_ar(2)/(t_ar(1)+t_ar(2));
	beta_crit=xthr*(t_ar(3)/t3c-alpha)
	
	for beta=0:5:50 % [0 1 10 50] % 0:0.1:10
	  tend=300;
	  r1=sim_dyn_nl(w,t_ar,tend,alpha,beta,xthr);
	  n2=length(r1)/2;
	  osc_ampl=mean(std(r1(:,n2:end)'))
	  res=[res [beta osc_ampl]'];
	  figure(2)
	  plot(res(1,:),res(2,:));
	  hold on
	  title('osc_ampl vs beta ')
	end  

case{31} % non-linear & damped osci.
	t_ar=[1 5 30]
	% (assuming gain=1)
% get_wcritimag([ 1 5 30]) = 0.237
% get_wcrit([ 1 5 30]) =0.89183
	  
	res=[];
	alpha=1;
	xthr=0.1;
	w=0.25;
	
	t3c=t_ar(1)*t_ar(2)/(t_ar(1)+t_ar(2)); % for w=0, I guess
	beta_crit=xthr*(t_ar(3)/t3c-alpha)
	
	%for beta=[0 1] % [0 1 10 50] % 0:0.1:10
	for beta =[0 1 2 3 4]
	  tend=300;
	  r1=sim_dyn_nl(w,t_ar,tend,alpha,beta,xthr);
	end  
	
case{4} % effective time-constants
	% network with recurrence, examine effective time-constnat as
	% one varies t3 from t3c*(1+eps) to t3i*10
	glob_functionplot=1
	w=0.99;
	tau1=1;
	t3c= get_tkcrit(w,[1 5])
	if (t3c <0) t3c=tau1; end
		
	t3ci= get_tkcritimag(w,[1 5])
   
	res=[];
	% need to make runtime adaptive
	%for t3=logspace(log10(t3c*1.05),log10(t3ci*10),50)
	for t3=20000
		[r1 taueff taueff_guess tauenv]=sim_dyn(w,[1 5 t3],1000);
		res=[res [t3 taueff taueff_guess tauenv]'];
	end
	
	taugoal=tau1/(1-w);
	res(5,:)=0*res(1,:)+taugoal;
	
	res(2:5,:) /= taugoal; % normalize, optional
	
	semilogx(res(1,:),res(2,:),';tau_eff normalized;');
	hold on
    plot(res(1,:),res(3,:),';tau-guess;');
	hold on
	plot(res(1,:),res(4,:),'r;tau-env;');
    plot(res(1,:),res(5,:),'g;tau-goal;');

	rest=res';
	save teff.dat rest

case{41} % effective time-constants, illustration.
	% network with recurrence, examine effective time-constnat as
	% one varies t3 from t3c*(1+eps) to t3i*10
	glob_functionplot=1
	%w=0.9;
	%t3c= get_tkcrit(w,[10 50])
	%r1=sim_dyn_demo(w,[10],1000);
	%r1=sim_dyn_demo(w,[10 50 t3c*1.5],1000);
	
	w=0.99;
	t3c= get_tkcrit(w,[10 50]);
	t3ci= get_tkcritimag(w,[10 50]);
	r1=sim_dyn_demo(w,[10],10000);
	r1=sim_dyn_demo(w,[10 50 t3c*1.5],10000);
	r1=sim_dyn_demo(w,[10 50 t3ci],10000);
	
	t3c*1.5
	t3ci
	
	error('stop')
	
	w=0.9;
	tau1=1;
	
	if (t3c <0) t3c=tau1; end
		
	t3ci= get_tkcritimag(w,[1 5])
   
	% need to make runtime adaptive
	%for t3=logspace(log10(t3c*1.05),log10(t3ci*10),50)
	for t3=200
		r1=sim_dyn_demo(w,[1 5 t3],100);
	end

case{5} % examine phase delays
	glob_functionplot=1
	
	t_ar=[1 5  10]
	tend=100
	w=0.5
	r1=sim_dyn(w,t_ar,tend);
	
	
	t_ar=[1 5 5 5 5 5 5 10]
	tend=100
	w=0.5
	r1=sim_dyn(w,t_ar,tend);
	
	
	t_ar=[1 5  100]
	tend=1000
	w=0.5
	r1=sim_dyn(w,t_ar,tend);
	

	
	
	
endswitch
