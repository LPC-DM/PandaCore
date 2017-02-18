#!/usr/bin/env python

import ROOT as root
import numpy as np
import root_numpy as rnp
from array import array
from PandaCore.Tools.Misc import *
from PandaCore.Tools.Load import Load
from os import getenv

Load('Drawers','HistogramDrawer')

class Process():
    def __init__(self,name,pt,custom_color=root.nProcesses):
        self.name = name
        self.process_type = pt
        if (custom_color==root.nProcesses):
            self.color = pt
        else:
            self.color = custom_color
        # public config - defaults
        self.dashed = False
        self.use_common_weight = True 
        self.use_common_cut = True
        self.additional_cut = '1==1'
        self.additional_weight = '1'
        self.tree_name = 'events'
        # private stuff
        self.__files = []
        return
    def add_file(self,fpath):
        self.__files.append(fpath)
    def read(self,variables,weights,cut):
        branches = variables + weights
        xarr = rnp.root2array(filenames=self.__files,
                              treename=self.tree_name,
                              branches=branches,
                              selection=cut)
        return xarr

def convert_name(n):
    rn = str(n)
    rn = rn.replace("/","Over");    
    rn = rn.replace("*","Times");   
    rn = rn.replace("+","Plus");    
    rn = rn.replace("-","Minus");   
    rn = rn.replace("TMath::","");  
    rn = rn.replace(")","");        
    rn = rn.replace("(","");        
    rn = rn.replace(".","_");       
    return rn


def divide_bin_width(h):
    nb = h.GetNbinsX()
    for ib in xrange(1,nb+1):
        val = h.GetBinContent(ib)
        err = h.GetBinError(ib)
        width = h.GetBinWidth(ib)
        h.SetBinContent(ib,val/width)
        h.SetBinError(ib,err/width)


def fix_underflow(h):
    h.SetBinContent(1,h.GetBinContent(0)+h.GetBinContent(1))

def fix_overflow(h):
    nbins = h.GetNbinsX()
    h.SetBinContent(nbins,h.GetBinContent(nbins+1)+h.GetBinContent(nbins))


class Distribution():
    def __init__(self,name,bin_list,is_variable,xlabel,ylabel,filename=None,ybounds=None):
        self.name = name
        if name=='1':
            # override
            is_variable = False
            bin_list = (0,2,1)
        self.is_variable = is_variable
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.filename = convert_name(filename if filename else name)
        self.ybounds = ybounds
        self.histograms = {} # so we can access all hists of this distribution
        self.systs = {}
        if is_variable:
            self.bins = array('f',bin_list)
            self.nbins = len(bin_list)-1
            self.hbase = root.TH1D('hbase_%s'%(self.filename),
                                   '',self.nbins,self.bins)
        else:
            self.xbounds = (bin_list[0],bin_list[1])
            self.nbins = bin_list[2]
            self.hbase = root.TH1D('hbase_%s'%(self.filename),
                                   '',self.nbins,
                                   self.xbounds[0],self.xbounds[1])
    def generate_hist(self,label):
        self.histograms[label] = self.hbase.Clone('h_%s_%s'%(self.filename,label))
        self.histograms[label].Sumw2()
        return self.histograms[label]
    def generate_syst(self,label):
        hup = self.hbase.Clone('h_%s_%s_Up'%(self.filename,label))
        hdown = self.hbase.Clone('h_%s_%s_Down'%(self.filename,label))
        hup.Sumw2(); hdown.Sumw2()
        self.systs[label] = (hup,hdown)
        return (hup,hdown)

def FDistribution(name,lo,hi,nbins,xlabel,ylabel,filename=None,ybounds=None):
    return Distribution(name=name,bin_list=(lo,hi,nbins),is_variable=False,
                        xlabel=xlabel,ylabel=ylabel,filename=filename,ybounds=ybounds)

def VDistribution(name,bins,xlabel,ylabel,filename=None,ybounds=None):
    return Distribution(name=name,bin_list=bins,is_variable=True,
                        xlabel=xlabel,ylabel=ylabel,filename=filename,ybounds=ybounds)

class Systematic():
    def __init__(self,name,weight_up,weight_down,color=1):
        self.name = name
        self.weight_up = weight_up
        self.weight_down = weight_down
        self.color = color
        self.hists = None
    def generate_weight(self,central,is_up):
        replacement = self.weight_up if is_up else self.weight_down
        if type(replacement) == str:
            return tTIMES(central,replacement)
        else:
            shift = str(central)
            for k,v in replacement.iteritems():
                shift.replace(k,v)
            return shift

class PlotUtility():
    def __init__(self):
        self.canvas = root.HistogramDrawer()
        members = [x for x in dir(self.canvas) if ('__' not in x)]
        for m in members:
            if callable(getattr(self.canvas,m)):
                self.__getPlotMethod(m)
            else:
                self.__getPlotMember(m)
        self.__processes = [] 
        self.__distributions = []
        self.__systematics = []
        # public configuration - defaults
        self.plot_label = '#it{CMS Preliminary}'
        self.do_underflow = True
        self.do_overflow = True
        self.signal_scale = 1
        self.cut = None
        self.mc_weight = None
        self.eventmod =  0
        self.eventnumber = 'eventNumber'
        return
    def __getPlotMethod(self,x):
        method = getattr(self.canvas,x)
        setattr(self,x, lambda *args : method(*args))
    def __getPlotMember(self,x):
        member = getattr(self.canvas,x)
        setattr(self,x, member)
    def add_process(self,p):
        self.__processes.append(p)
    def add_distribution(self,d):
        self.__distributions.append(d)
    def add_systematic(self,name,up,down,color=1):
        self.__systematics.append(Systematic(name,up,down,color))
    def Draw(self):
        # override the thing that's bound in __init__
        PError('plot_utility.PlotUtility.Draw','Not implemented!')
    def draw_all(self,outdir):
        if not self.canvas.HasLegend():
            self.canvas.InitLegend()

        f_out = root.TFile(outdir+'hists.root','UPDATE')
        if f_out.IsZombie():
            f_out.close()
            f_out = root.TFile(outdir+'hists.root','RECREATE')
        f_buffer = root.TFile('/tmp/%s/buffer_%i.root'%(getenv('USER'),root.gSystem.GetPid()),'RECREATE')
        f_buffer.cd()

        variables = []
        for dist in self.__distributions:
            for syst in self.__systematics:
                syst.hists = dist.generate_syst(syst.name)

            for proc in self.__processes:
                dist.generate_hist(proc.name)

            variables.append(dist.name)

        # loop through each process
        for proc in self.__processes:
            # figure out the nominal weight and cut strings
            final_weight = '1'
            final_cut = self.cut
            if (proc.process_type!=root.kData and proc.use_common_weight):
                final_weight = self.mc_weight
            if (proc.process_type<=root.kSignal3 and proc.process_type!=root.kData):
                final_weight = tTIMES(final_weight,str(self.signal_scale))
            if self.eventmod:
                if (p.process_type==root.kData):
                    final_cut = tAND(final_cut,'(%s%%%i)==0'%(self.eventnumber,self.eventmod))
                else:
                    final_weight = tTIMES(final_weight,str(1./self.eventmod))
            final_cut = tAND(final_cut,proc.additional_cut)
            final_weight = tTIMES(final_weight,proc.additional_weight)

            weight_map = {'nominal' : final_weight}
            weights = [final_weight]
            if proc.process_type!=root.kData:
                for syst in self.__systematics:
                    if proc.use_common_weight:
                        up_weight = syst.generate_weight(final_weight,True)
                        down_weight = syst.generate_weight(final_weight,False)
                    else:
                        continue
                    if up_weight in weight_map or down_weight in weight_map:
                        continue # do not request duplicate branches
                    weight_map['%s_Up'%(syst.name)] = up_weight
                    weight_map['%s_Down'%(syst.name)] = down_weight
                    weights.append(up_weight)
                    weights.append(down_weight)

            xarr = proc.read(variables,weights,final_cut)

            for dist in self.__distributions:
                vals = xarr[dist.name]
                weights_nominal = xarr[weight_map['nominal']]
                rnp.fill_hist(hist=dist.histograms[proc.name],
                              array=vals,
                              weights=weights_nominal)
                if proc.process_type!=root.kData:
                    for syst in self.__systematics:
                        if proc.use_common_weight:
                            weights_up = xarr[weight_map['%s_Up'%syst.name]]
                            weights_down = xarr[weight_map['%s_Down'%syst.name]]
                        else:
                            weights_up = weights_nominal
                            weights_down = weights_nominal
                        rnp.fill_hist(hist=dist.systs[syst.name][0],
                                      array=vals,
                                      weights=weights_up)
                        rnp.fill_hist(hist=dist.systs[syst.name][1],
                                      array=vals,
                                      weights=weights_down)

        # everything is filled, now draw the histograms!
        for dist in self.__distributions:
            if dist.name=="1":
                totals = {'bg':0, 'sig':0, 'data':0}
                errs = {'bg':0, 'sig':0, 'data':0}
                table_format = '%-25s | %15f +/- %15f'

                table = ['%-25s | %15s +/- %15s'%('Process','Yield','Stat. unc.')]
                table.append("=================================================================")

                for proc in self.__processes:
                    h = dist.histograms[proc.name]
                    integral = h.GetBinContent(1)
                    error = h.GetBinError(1)
                    table.append(table_format%(proc.name,integral,error))
                    if proc.process_type==root.kData:
                        proc_label = 'data'
                    elif proc.process_type<=root.kSignal3:
                        proc_label = 'sig'
                    else:
                        proc_label = 'bg'
                    totals[proc_label] += integral
                    errs[proc_label] += pow(error,2)

                for k,v in errs.iteritems():
                    errs[k] = np.sqrt(v)

                table.append("=================================================================")
                table.append(table_format%('MC(bkg)',totals['bg'],errs['bg']))
                table.append(table_format%('MC(sig)',totals['sig'],errs['sig']))
                table.append(table_format%('Data',totals['data'],errs['data']))

                table.append("=================================================================")
                if totals['bg']:
                    table.append('S/B=%.3f, S/sqrtB=%.3f'%(totals['sig']/totals['bg'],
                                                           totals['sig']/np.sqrt(totals['bg'])))

                with open(outdir+'yields.txt','w') as fyields:
                    fyields.write('\n'.join(table))

                for t in table:
                    PInfo('plot_utility.PlotUtility.Dump',t)

            self.canvas.Reset(False)
            for proc in self.__processes:
                h = dist.histograms[proc.name]
                if self.do_overflow:
                    fix_overflow(h)
                if self.do_underflow:
                    fix_underflow(h)
                if dist.is_variable:
                    divide_bin_width(h)
                if dist.ybounds:
                    h.SetMinimum(ybounds[0])
                    h.SetMaximum(ybounds[1])
                h.GetXaxis().SetTitle(dist.xlabel)
                h.GetYaxis().SetTitle(dist.ylabel)
                h.SetTitle('')
                if self.canvas.IsStack():
                    h.SetLineWidth(2)
                if (proc.process_type<=root.kSignal3 
                       and proc.process_type!=root.kData 
                       and self.signal_scale!=1):
                    self.canvas.AddHistogram(h,'%.1f#times%s'%(self.signal_scale,proc.name),proc.process_type)
                else:
                    self.canvas.AddHistogram(h,proc.name,proc.process_type)
            for syst in self.__systematics:
                hup,hdown = dist.systs[syst.name]
                for h in [hup,hdown]:
                    if self.do_overflow:
                        fix_overflow(h)
                    if self.do_underflow:
                        fix_underflow(h)
                    if dist.is_variable:
                        divide_bin_width(h)
                    h.SetLineWidth(3)
                    h.SetLineColor(syst.color)
                self.canvas.AddSystematic(hup,'hist',syst.name)
                self.canvas.AddSystematic(hdown,'hist')

            # output the canvas
            self.canvas.Logy(False)
            self.canvas.Draw(outdir,dist.filename)
            self.canvas.ClearLegend()
            self.canvas.Logy(True)
            self.canvas.Draw(outdir,dist.filename+'_logy')







