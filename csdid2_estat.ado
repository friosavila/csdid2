program clrreturn, rclass
        exit
end

program adde, eclass
        ereturn `0'
end

program addr, rclass
		return add
        return `0'
end

program adds, sclass
        sreturn `0'
end

 program csdid2_estat, sortpreserve  
	version 14
		syntax anything, [*]
        capture mata:csdid
		if _rc!=0   error 301
		
		gettoken key rest : 0, parse(", ")
		
		if inlist("`key'","attgt","simple","pretrend","group","calendar","event","cevent") {
			csdid_do `key' `rest'
			
		}
		else {
		    display in red "Option `key' not recognized"
			error 199
		}
		
end
 program csdid_do, rclass
	syntax namelist(min=1 max=1 ), [ estore(name) ///
					esave(name) replace ///
					post level(int 95)  ///
					WBOOT				///
					WBOOT1(str)			///
					reps(integer 999) 	///
					rseed(string) 		///
					wbtype(string)		///
					rgroup(numlist)     ///
					rcalendar(numlist) /// 
					revent(numlist)    ///
					max_mem(real 1)  plot  * ]
	
	// confirm csdid exists and if csdidstat=csdid_estat()
	local key `namelist'
	capture mata:csdidstat
	if _rc!=0 mata:csdidstat=csdid_estat()
	
	// check Keys
	
	if "`key'"=="pretrend" local ktype = 0
	if "`key'"=="attgt"    local ktype = 1
	if "`key'"=="simple"   local ktype = 2
	if "`key'"=="group"    local ktype = 3
	if "`key'"=="calendar" local ktype = 4
	if "`key'"=="event"    local ktype = 5 
	if "`key'"=="cevent"   local ktype = 6
	
	/// initialize
	if "`rseed'"!="" set seed `rseed'
	
	mata: csdidstat.cilevel = `level'/100
	mata: csdidstat.bwtype  = 1      
	mata: csdidstat.reps    = `reps'
	mata: csdidstat.max_mem = `max_mem'
	mata: csdidstat.range.selgvar = J(0,0,.)
	mata: csdidstat.range.seltvar = J(0,0,.)
	mata: csdidstat.range.selevent= J(0,0,.)
	
	/// Check if we have to make sample selection
	
	if "`rgroup'"!="" {
		 numlist "`rgroup'", int
		 mata:csdidstat.range.selgvar=csdidstat.rtokens("`r(numlist)'")
	}
	if "`rcalendar'"!="" {
		 numlist "`rcalendar'", int
		 mata:csdidstat.range.seltvar=csdidstat.rtokens("`r(numlist)'")
	}
	if "`revent'"!="" {
		 numlist "`revent'", int
		 mata:csdidstat.range.selevent=csdidstat.rtokens("`r(numlist)'")
	}
	
	if `ktype'>0 		mata: csdidstat.test_type  = `ktype'      
	else {
		mata: csdidstat.pretrend(csdid)
		display "Pre-trend test"
		display "H0: All ATTGT=g for all T<G"
		display "chi2(`r(df)') = " %10.4f scalar(chi2_)
		display "p-value  = " %10.4f scalar(pchi2_)	
		return   scalar chi2  = scalar(chi2_)
		return   scalar pchi2 = scalar(pchi2_)	
		return   scalar df    = scalar(df_)	
		exit
	}
	
	if "`wboot'`wboot1'"=="" {
		mata:csdidstat.atts_asym(csdid)
		
		capture:est store `lastreg'	
		ereturn clear
		return matrix b = _bb, copy
		return matrix V = _vv, copy
		
		adde post _bb _vv
		adde local cmd 	   csdid2
		adde local estat_cmd csdid2_estat	
		//syntax namelist(min=1 max=1 ), [*]
		adde local cmdline estat  `0'
		adde local agg     `key'
		adde local aggt     `ktype'
 		
		if "`estore'"!="" est store `estore'
		if "`esave'" !="" est save  `esave', `replace'
		_coef_table, level(`level')
		matrix rtb=r(table)
		
		if "`post'"=="" qui:est restore `lastreg'

		return matrix table = rtb, copy
		return local agg  `key'
		
	}
	else {
		mata:csdidstat.atts_wboot(csdid)
		
		capture:est store `lastreg'	
		ereturn clear
		return matrix table = _table, copy
 		tempname bb
		matrix `bb' = _table[1,....]
		return matrix b = `bb', copy
		adde post `bb'
		adde matrix table1 = _table, copy
		adde local vcetype WBoot
		adde local cmd 	   csdid2
		adde local estat_cmd csdid2_estat
		
		//syntax namelist(min=1 max=1 ), [*]
		adde local cmdline estat `0'
		adde local agg     `key'
		
		csdid_tablex, `diopts' level(`level')	
		display "WildBootstrap Standard errors"	_n ///
				"with `reps' Repetitions"
		matrix rtb = r(table)
		
 
	
		if "`post'"=="" capture:qui:est restore `lastreg'
		return matrix table = rtb, copy
		return local agg  `key'
		
	}

		
	if "`plot'"!="" {	    
	    csdid_plots ,  ktype(`ktype') `options'
	}
	capture matrix drop rtb 
end

program csdid_plots
	syntax, [style(string) * ktype(int 5)] 	
	if "`style'"=="" local style rspike
	if !inlist("`style'","rspike","rarea","rcap","rbar") {
	    display as erro "Style `style' not allowed"
	    error 1
	}
	if `ktype'==5 {
		tempvar t b ll uu
		mata:event_p("`t' `b' `ll' `uu'")
		
		csdid_plot_eventx 	`t' `b' `ll' `uu', style(`style') `options'		
	}
	else if `ktype'==3 | `ktype'==4 {
	    // Group Calendar
		tempvar t b ll uu
		mata:other_p("`t' `b' `ll' `uu'")
		
		csdid_plot_other `t' `b' `ll' `uu', style(`style') `options' ktype(`ktype')		
	}
	else {
	    display in red "Plot option not allowed"
	}
end

		



program csdid_plot_eventx 
	syntax varlist, style(str) [pstyle1(string asis) pstyle2(string asis) ///
								 color1(string asis)  color2(string asis) ///
								 LWidth(string asis) barwidth(string asis) ///
								 xtitle(passthru)     ytitle(passthru) ///
								 legend(passthru)  * ]
	gettoken t rest:varlist
	gettoken b rest:rest
	gettoken ll rest:rest 
	gettoken uu rest:rest 
	** defaults
	if "`xtitle'"=="" local xtitle xtitle("Periods to Treatment")
	if "`ytitle'"=="" local ytitle ytitle("ATT")
	if "`pstyle1'"=="" local pstyle1 p1
	if "`pstyle2'"=="" local pstyle2 p2
	if `"`color1'"'=="" local color1 %40
	if `"`color2'"'=="" local color2 %40
	if "`legend'"=="" local legend legend(order(1 "Pre-treatment" 3 "Post-treatment"))	
	
	if "`style'"=="rspike" {
		if "`lwidth'"    =="" local lwidth 3		
		two   `style' `ll' `uu' `t'   if `t'<0 , pstyle(`pstyle1') color(`color1') lw(`lwidth') || ///
			  scatter  `b'      `t'   if `t'<0 , pstyle(`pstyle1') || ///
			  `style' `ll' `uu' `t'   if `t'>=0, pstyle(`pstyle2') color(`color2') lw(`lwidth') || ///
			  scatter  `b'      `t'   if `t'>=0, pstyle(`pstyle2') , ///
			  `legend'  `xtitle' `ytitle' ///
			  yline(0 , lp(dash) lcolor(black))   `options'	  
	}	  
	
	else if "`style'"=="rarea" {
	 	 if "`lwidth'"    =="" local lwidth 0
	two   (`style'  `ll' `uu' `t'   if `t'<0 , pstyle(`pstyle1') color(`color1') lw(`lwidth') ) || ///
		  (scatter  `b'     `t'   if `t'<0   , pstyle(`pstyle1') connect(l) ) || ///
		  (`style'  `ll' `uu' `t'   if `t'>=0, pstyle(`pstyle1') color(`color2') lw(`lwidth') ) || ///
		  (scatter  `b'     `t'   if `t'>=0  , pstyle(`pstyle2') connect(l) ), ///
		  `legend'  `xtitle' `ytitle' ///
		  yline(0 , lp(dash) lcolor(black))   `options'
	}
	
	else if "`style'"=="rcap" {
	    if "`lwidth'"    =="" local lwidth 1
	two   (`style'  `ll' `uu' `t'   if `t'<0, pstyle(`pstyle1') color(`color1') lw(`lwidth') ) || ///
		  (scatter  `b'      `t'   if `t'<0 , pstyle(`pstyle1') connect(l) ) || ///
		  (`style' `ll' `uu' `t'   if `t'>=0, pstyle(`pstyle2') color(`color2') lw(`lwidth') ) || ///
		  (scatter  `b'      `t'   if `t'>=0, pstyle(`pstyle2') connect(l) ), ///
		   `legend'  `xtitle' `ytitle' ///
		  yline(0 , lp(dash) lcolor(black))    `options'
	}
	
	else if "`style'"=="rbar" {
	    if "`lwidth'"    =="" local lwidth 0
		if "`barwidth'"  =="" local barwidth 0.5
	two   (`style'  `ll' `uu' `t'   if `t'<0 , pstyle(`pstyle1') color(`color1') lw(`lwidth') barwidth(`barwidth') ) || ///
		  (scatter  `b'      `t'   if `t'<0  , pstyle(`pstyle1') connect(l) ) || ///
		  (`style' `ll' `uu' `t'   if `t'>=0 , pstyle(`pstyle2') color(`color1') lw(`lwidth') barwidth(`barwidth') ) || ///
		  (scatter  `b'      `t'   if `t'>=0 , pstyle(`pstyle2') connect(l) ), ///
		  `legend'   `xtitle' `ytitle' ///
		  yline(0 , lp(dash) lcolor(black))  `options'
	} 
 
end

program csdid_plot_other
	syntax varlist, style(str) [ktype(int 3) * ///
								 pstyle(passthru) color(passthru)  ///
								 LWidth(passthru) barwidth(passthru) ///
								 xtitle(passthru) ytitle(passthru) ///
								 legend(passthru)]
	gettoken t rest:varlist
	gettoken b rest:rest
	gettoken ll rest:rest 
	gettoken uu rest:rest 
	
 
	if `"`pstyle'"'=="" local pstyle pstyle(p1)
	if `"`color'"'=="" local color color(%40)
	
	tempvar tt 
	encode `t', gen(`tt')
	
	if `"`xtitle'"'=="" & `ktype' ==3 local xtitle xtitle("Group")
	else if `"`xtitle'"'=="" & `ktype' ==4 local xtitle xtitle("Calendar")
	if `"`ytitle'"'=="" local ytitle ytitle("ATT")
	

	if "`style'"=="rspike" {
	    if "`lwidth'"    =="" local lwidth lwidth(3)
	two   (`style'  `ll' `uu' `tt'   ,  `pstyle'  `color' `lwidth' ) || ///
		  (scatter  `b'      `tt'    ,  `pstyle'    ) , ///
		  legend(off) `xtitle'  `ytitle'  ///
		  yline(0 , lp(dash) lcolor(black)) `title' `name' xlabel(, val) `options'
	}	  
	
	if "`style'"=="rarea" {
	    if "`lwidth'"    =="" local lwidth lwidth(0)
	two   (`style'  `ll' `uu' `tt'   , `pstyle'  `color' `lwidth' ) || ///
		  (scatter  `b'     `tt'     , `pstyle'  ) , ///
		  legend(off) `xtitle'  `ytitle'  ///
		  yline(0 , lp(dash) lcolor(black)) `title' `name' xlabel(, val) `options'
	}
	
	if "`style'"=="rcap" {
	    if "`lwidth'"    =="" local lwidth lwidth(1)
	two   (`style'  `ll' `uu' `tt'   , `pstyle'  `color' `lwidth' ) || ///
		  (scatter  `b'     `tt'   	 , `pstyle'  ) , ///
		  legend(off) `xtitle'  `ytitle'  ///
		  yline(0 , lp(dash) lcolor(black)) `title' `name' xlabel(, val) `options'
	}
	
	if "`style'"=="rbar" {
	    if "`lwidth'"    =="" local lwidth lwidth(0)
		if "`barwidth'"  =="" local barwidth 0.5
	two   (`style'  `ll' `uu' `tt'   	 , `pstyle'  `color' `lwidth') || ///
		  (scatter  `b'     `tt'   	 , `pstyle'  ) , ///
		  legend(off) `xtitle'  `ytitle'  ///
		  yline(0 , lp(dash) lcolor(black)) `title' `name' xlabel(, val) `options'
	}
end

mata:
 	void event_p( string scalar newvars){
	    real   matrix tbl, ntbl2
		string matrix ntbl
	    tbl = st_matrix("rtb")		
		ntbl = st_matrixcolstripe("rtb")
		ntbl = usubinstr(ntbl,"tp","+",.)
		ntbl = usubinstr(ntbl,"tm","-",.)	
		ntbl2= strtoreal(ntbl)	
		tbl  = tbl[(1,5,6),]'	
		tbl  = select(tbl,(ntbl2[,2]:!=.))		
		ntbl2= select(ntbl2[,2],(ntbl2[,2]:!=.))
        real matrix ss
 		ss= _st_addvar("double",tokens(newvars))
 		st_store((1::rows(tbl)) ,tokens(newvars),(ntbl2,tbl))	
	}
 
	void other_p(string scalar newvars){
	    real   matrix tbl
		string matrix ntbl
	    tbl  = st_matrix("rtb")		
		ntbl = st_matrixcolstripe("rtb")
		ntbl = ntbl [,2]
		tbl  = tbl[(1,5,6),]'	
		string matrix tnv
		tnv = tokens(newvars)
		real matrix ss
		ss= _st_addvar(sprintf("str%f",max(strlen(ntbl))),tnv[1])
		ss= _st_addvar("double",tnv[2..4])
		st_sstore((1::rows(tbl)) ,tnv[1],ntbl)	
		st_store((1::rows(tbl)) ,tnv[2..4],tbl)	
	}
end
 

program csdid_plot_event 
	syntax varlist, style(str) [title(passthru) name(passthru) ///
								ytitle(passthru) xtitle(passthru)	///
								legend(str) * ]
	gettoken t rest:varlist
	gettoken b rest:rest
	gettoken ll rest:rest 
	gettoken uu rest:rest 
	** defaults
	if "`xtitle'"=="" local xtitle xtitle("Periods to Treatment")
	if "`ytitle'"=="" local ytitle ytitle("ATT")
	
		
	if "`style'"=="rspike" {
	two   rspike  `ll' `uu' `t'   if `t'<0 , pstyle(p1) color(%40) lw(3) || ///
		  scatter  `b'      `t'   if `t'<0 , pstyle(p1) || ///
		  rspike  `ll' `uu' `t'   if `t'>=0, color(%40) pstyle(p2) lw(3) || ///
		  scatter  `b'      `t'   if `t'>=0, pstyle(p2) , ///
		  legend(order(1 "Pre-treatment" 3 "Post-treatment") `legend' ) `xtitle' `ytitle' ///
		  yline(0 , lp(dash) lcolor(black)) `title' `name'  `options'
	}	  
	
	if "`style'"=="rarea" {
	two   (rarea  `ll' `uu' `t'   if `t'<0 , pstyle(p1) color(%40) lw(0) ) || ///
		  (scatter  `b'     `t'   if `t'<0 , pstyle(p1) connect(l) ) || ///
		  (rarea  `ll' `uu' `t'   if `t'>=0, color(%40) pstyle(p2) lw(0) ) || ///
		  (scatter  `b'     `t'   if `t'>=0, pstyle(p2) connect(l) ), ///
		  legend(order(1 "Pre-treatment" 3 "Post-treatment") `legend' ) `xtitle' `ytitle' ///
		  yline(0 , lp(dash) lcolor(black))  `title' `name'  `options'
	}
	
	if "`style'"=="rcap" {
	two   (rcap  `ll' `uu' `t'   if `t'<0, pstyle(p1) color(%60) lw(1) ) || ///
		  (scatter  `b'      `t'   if `t'<0 , pstyle(p1) connect(l) ) || ///
		  (rcap `ll' `uu' `t'   if `t'>=0, color(%60) pstyle(p2) lw(1) ) || ///
		  (scatter  `b'      `t'   if `t'>=0, pstyle(p2) connect(l) ), ///
		  legend(order(1 "Pre-treatment" 3 "Post-treatment")  `legend') `xtitle' `ytitle' ///
		  yline(0 , lp(dash) lcolor(black))  `title' `name' `options'
	}
	
	if "`style'"=="rbar" {
	two   (rbar  `ll' `uu' `t'   if `t'<0, pstyle(p1) color(%60) lw(0) barwidth(0.5) ) || ///
		  (scatter  `b'      `t'   if `t'<0 , pstyle(p1) connect(l) ) || ///
		  (rbar `ll' `uu' `t'   if `t'>=0, color(%60) pstyle(p2) lw(0) barwidth(0.5) ) || ///
		  (scatter  `b'      `t'   if `t'>=0, pstyle(p2) connect(l) ), ///
		  legend(order(1 "Pre-treatment" 3 "Post-treatment") `legend' ) `xtitle' `ytitle' ///
		  yline(0 , lp(dash) lcolor(black))  `title' `name' `options'
	} 
end



program csdid_plot_group
	syntax varlist, style(str) [title(passthru) name(passthru)	///
								ytitle(passthru) xtitle(passthru) * ]
	gettoken t rest:varlist
	gettoken b rest:rest
	gettoken ll rest:rest 
	gettoken uu rest:rest 
		
	qui:levelsof `t', local(tlev)
	local tlb: value label `t'
	local xlab 0 " "
	foreach i of local tlev {
	    local j = `j'+1
	    local xlab `xlab' `i' "`:label `tlb' `i''"
	}
	
	if "`xtitle'"=="" local xtitle xtitle("Groups")
	if "`ytitle'"=="" local ytitle ytitle("ATT")
	
	local xlab `xlab' `=`j'+1' " "
	
	if "`style'"=="rspike" {
	two   (`style'  `ll' `uu' `t'   , pstyle(p1) color(%40) lw(3) ) || ///
		  (scatter  `b'      `t'   , pstyle(p1)    ) , ///
		  legend(off) `xtitle'  `ytitle'  ///
		  yline(0 , lp(dash) lcolor(black)) `title' `name' xlabel(`xlab') `options'
	}	  
	
	if "`style'"=="rarea" {
	two   (`style'  `ll' `uu' `t'   	 , pstyle(p1) color(%40) lw(0) ) || ///
		  (scatter  `b'     `t'   	 , pstyle(p1)  ) , ///
		  legend(off) `xtitle'  `ytitle'  ///
		  yline(0 , lp(dash) lcolor(black)) `title' `name' xlabel(`xlab') `options'
	}
	
	if "`style'"=="rcap" {
	two   (`style'  `ll' `uu' `t'   	 , pstyle(p1) color(%40)  lw(1) ) || ///
		  (scatter  `b'     `t'   	 , pstyle(p1)  ) , ///
		  legend(off) `xtitle'  `ytitle'  ///
		  yline(0 , lp(dash) lcolor(black)) `title' `name' xlabel(`xlab') `options'
	}
	
	if "`style'"=="rbar" {
	two   (`style'  `ll' `uu' `t'   	 , pstyle(p1) color(%40) lw(0) barwidth(0.5) ) || ///
		  (scatter  `b'     `t'   	 , pstyle(p1)  ) , ///
		  legend(off) `xtitle'  `ytitle'  ///
		  yline(0 , lp(dash) lcolor(black)) `title' `name' xlabel(`xlab') `options'
	}
end


program csdid_plot_calendar
	syntax varlist, style(str) [title(passthru) name(passthru)	///
								ytitle(passthru) xtitle(passthru) * ]
	gettoken t rest:varlist
	gettoken b rest:rest
	gettoken ll rest:rest 
	gettoken uu rest:rest 
		
	qui:levelsof `t', local(tlev)
	local tlb: value label `t'
	local xlab 0 " "
	foreach i of local tlev {
	    local j = `j'+1
	    local xlab `xlab' `i' "`:label `tlb' `i''"
	}
	
	if "`xtitle'"=="" local xtitle xtitle("Periods")
	if "`ytitle'"=="" local ytitle ytitle("ATT")
	

	local xlab `xlab' `=`j'+1' " "
	
	if "`style'"=="rspike" {
	two   (`style'  `ll' `uu' `t'   , pstyle(p1) color(%40) lw(3) ) || ///
		  (scatter  `b'      `t'   , pstyle(p1)    ) , ///
		  legend(off) `xtitle'  `ytitle'  ///
		  yline(0 , lp(dash) lcolor(black)) `title' `name' xlabel(`xlab') `options'
	}	  
	
	if "`style'"=="rarea" {
	two   (`style'  `ll' `uu' `t'   	 , pstyle(p1) color(%40) lw(0) ) || ///
		  (scatter  `b'     `t'   	 , pstyle(p1)  ) , ///
		  legend(off) `xtitle'  `ytitle'  ///
		  yline(0 , lp(dash) lcolor(black)) `title' `name' xlabel(`xlab') `options'
	}
	
	if "`style'"=="rcap" {
	two   (`style'  `ll' `uu' `t'   	 , pstyle(p1) color(%40)  lw(1) ) || ///
		  (scatter  `b'     `t'   	 , pstyle(p1)  ) , ///
		  legend(off) `xtitle'  `ytitle'  ///
		  yline(0 , lp(dash) lcolor(black)) `title' `name' xlabel(`xlab') `options'
	}
	
	if "`style'"=="rbar" {
	two   (`style'  `ll' `uu' `t'   	 , pstyle(p1) color(%40) lw(0) barwidth(0.5) ) || ///
		  (scatter  `b'     `t'   	 , pstyle(p1)  ) , ///
		  legend(off) `xtitle'  `ytitle'  ///
		  yline(0 , lp(dash) lcolor(black)) `title' `name' xlabel(`xlab') `options'
	}
end

program define tsvmat2, return
        syntax anything, name(string)
        version 7
		 
        local nx = rowsof(matrix(`anything'))
        local nc = colsof(matrix(`anything'))
        ***************************************
        // here is where the safegards will be done.
        if _N<`nx' {
            display as result "Expanding observations to `nx'"
                set obs `nx'
        }
        // here we create all variables
        foreach i in `name' {
			local j = `j'+1
			qui:gen `type' `i'=matrix(`anything'[_n,`j'])			
        }
        // here is where they are renamed.

end


program csdid_tablex, rclass 
	syntax [, level(int `c(level)') noci cformat(string) sformat(string) *]

	_get_diopts diopts rest, `options'

	local cf %9.0g  
	local pf %5.3f
	local sf %7.2f

	if ("`cformat'"!="") {
			local cf `cformat'
	}
	if ("`sformat'"!="") {
			local sf `sformat'
	}
***hack to get max
 tempname tablex
 matrix `tablex' = _table'
 local namelist : colname `tablex'
 local wdt=0
 foreach i of local namelist {
 	if length("`i'")>`wdt' local wdt = length("`i'")+3
 }
 if `wdt'<15 local wdt = 12
***
        tempname mytab b se z t  ll ul cimat rtab
        .`mytab' = ._tab.new, col(6) lmargin(0)
        .`mytab'.width    `wdt'   |12    12     8         12    12
        .`mytab'.titlefmt  .     .     .   %6s       %24s     .
        .`mytab'.pad       .     2     1     0          3     3
        .`mytab'.numfmt    . %9.0g %9.0g %7.2f    %9.0g %9.0g

		
		local stat t 
		
        local namelist : rowname `tablex'		
        local eqlist : roweq `tablex'
        local k : word count `namelist'
		local knew = `k'
		matrix `rtab' = J(9, `k', .)
		matrix `cimat'= `tablex'
		* pvalue
		matrix rownames `rtab' = b se t p ll ul df crit eform
		matrix colnames `rtab' = `namelist'
		forvalues i = 1/`k' {
		    local kxc: word `i' of `eqlist'
			if ("`kxc'"=="wgt") {
				local knew = `knew' -1
			}
			matrix `rtab'[1,`i'] = `cimat'[`i',1]
			matrix `rtab'[2,`i'] = `cimat'[`i',2]
			matrix `rtab'[3,`i'] = `cimat'[`i',3]
			matrix `rtab'[5,`i'] = `cimat'[`i',4]
			matrix `rtab'[6,`i'] = `cimat'[`i',5]
			matrix `rtab'[8,`i'] = `cimat'[`i',6]
		}
        .`mytab'.sep, top
        if `:word count `e(depvar)'' == 1 {
                local depvar "`e(depvar)'"
        }
        .`mytab'.titles "`depvar'"                      /// 1
                        " Coefficient"                  /// 2
                        "Std. err."                     /// 3
                        "`stat'"                        /// 4   "P>|`stat'|"                    /// 5
                        "[`level'% conf. interval]" ""  //  6 7
		
        forvalues i = 1/`knew' {
                local name : word `i' of `namelist'
                local eq   : word `i' of `eqlist'
                if ("`eq'" != "_") {
                        if "`eq'" != "`eq0'" {
                                .`mytab'.sep
                                local eq0 `"`eq'"'
                                .`mytab'.strcolor result  .  .  .  .    .
                                .`mytab'.strfmt    %-12s  .  .  .  .    .
                                .`mytab'.row      "`eq'" "" "" "" ""  ""
                                .`mytab'.strcolor   text  .  .  .  .    .
                                .`mytab'.strfmt     %12s  .  .  .  .    .
                        }
                        local beq "[`eq']"
                }
                else if `i' == 1 {
                        local eq
                        .`mytab'.sep
                }
				
                scalar `b' = `cimat'[`i',1]
				scalar `se' = `cimat'[`i',2]
				scalar `t' = `cimat'[`i',3]
				scalar `ll'   = `cimat'[`i',4]
				scalar `ul'   = `cimat'[`i',5]
                .`mytab'.row    "`name'"                ///
                                `b'         ///
                                `se'        ///
                                `t'                     /// `p'  ///
                                `ll' `ul'
        }
        .`mytab'.sep, bottom
		return matrix table = `rtab'
end