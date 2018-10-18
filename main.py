#!/usr/bin/python
from __future__ import print_function
import os, time
from copy import deepcopy as dcp
def Type(x):
	return str(type(x)).split(" ")[1][1:-2]
def Str2Num(Str):
	try:
		return int(Str)
	except:
		try:
			return float(Str)
		except:
			return None

class Attibute:
	def __init__(self,Name="Omari",data=[]):
		self.name=Name
		self.value=[]
		self.__index__=[]
		self.__meta__=[]
		#  
		self.Import(data)
	def __repr__(self):
		return "OMARI"
		#self.Return()
	def Len(self):
		return len(self.value)
	
	def __Copy__(self):
		if len(self.__meta__)==0:self.__meta__=dcp(self.value)

	def Reset(self):
		if len(self.__meta__) !=00:
			self.value=dcp(self.__meta__)
			self.__meta__=[]   
	
	def Update(self):
		self.__meta__=dcp(self.value)
		self.Reset()
		return self

	def Import(self,Data):
		__Type__ = Type(Data)
		if __Type__=="list":
			self.value+=Data
		else:
			self.value.append(Data)
		return self

	def Return(self):
		return self.value
	
	# ******************
	def Slice(self,I=[]):#start=None,end=None,step=1):
		self.__Copy__()
		self.value=[self.value[i] for i in I] #self.value(I)
		return self

	def Filter(self,Function,Case=0):
		self.__Copy__()
		A = filter(Function,self.value)
		I=[]; initial=0
		for j in range(len(A)):
			for i in range(initial,self.Len()):
				if A[j]==self.value[i]:
					initial=i+1
					I.append(i)
					break
		self.value=A
		self.__index__=I
		return self

	def Map(self,Function):
		self.__Copy__()
		self.value = map(Function,self.value)
		return self
	def ExecutueFunction(self,Function):
		self.__Copy__()
		self.value=Function(dcp(self.value))
		return self
	def Convert2Numbers(self):
		self.Map(Str2Num)
		return self
		

# =================================================================================
class database:
	def __init__(self):
		self.name = "omari"
		self.createtime=0
		self.attributes=[]
		self.__selected__att__=[]
		#self.__meta__={}
		self.__export_to_globels__=False
		#
		self.__functions__={}
		self.__methode_att__=[]
		pass
	def __Copy__(self):
		#if len(self.__meta__)==0:self.__meta__=dcp(self.value)
		pass
	def ShortCut(self,On=True,Case=0):
		"""
			Export the class variables to the globals valible
			Case 0: Class
			Case 1: Str
		"""
		if Case == 0:
			for att in self.attributes:globals()[att]=self.__dict__[att]
		if Case == 1:
			for att in self.attributes:globals()[att]=self.__dict__[att]
		pass
	# ----------------------------------------------------------------------------------------
	
	def Reset(self):
		for att in self.attributes:self.__dict__[att].Reset()
		pass

	def Update(self):
		pass

	# ----------------------------------------------------------------------------------------
	def AddAtribute(self,Name="OMARI",Data=[]):
		if Name in self.__dict__.keys():
			pass	
		else:
			self.__dict__[Name]=Attibute(Name,Data)
			self.attributes.append(Name)
			pass
		pass

	def RestAttriute(self,Name):
		if Name in self.__dict__.keys():self.__dict__[Name]=Attibute(Name)
		pass

	def Select(self,*argu,**kwargs):
		if "Reset" in kwargs:
			if kwargs["Reset"]:
				self.__selected__att__=[]
		Case=False
		for att in argu:
			if att in self.attributes:
				Case=True
				self.__selected__att__=att
		if len(argu)==0 and not Case:
			self.__selected__att__=self.attributes
		pass

	def __use_sel_att__(self):
		if len(self.__selected__att__)==0:
			return self.attributes
		return self.__selected__att__
	
	def Del(self):
		pass
	def GenerateIndexArray(self):
		self.AddAtribute(Name="__index__",Data=list(range(self.__dict__[self.attributes[0]].Len())))
	# ----------------------------------------------------------------------------------------
	def Function(self,Function,MetaName=""):
		#self.__meta__[MetaName]=Function(self)
		return self

	# -----------------------------------------------------------------------------------------
	def Slice(self,I=[]):
		for att in self.__use_sel_att__():self.__dict__[att].Slice(I)
		return self

	def Filter(self,Condition,AttibuteName=""):
		I=self.__dict__[AttibuteName].Filter(Condition).__index__
		self.Slice(I)
		return self
	# ------------------------------------------------------------------------------------------
	def Return(self):
		return [self.__dict__[att].Return() for att in self.__use_sel_att__()]
	
	# ------------------------------------------------------------------------------------------
	#@staticmethod
	def Functions(self,Function,Name):
		#@classmethod
		class Fun1:
			def __init__(self,outer):
				self.outer=outer
			def Add(self):
				if Name in self.outer.__dict__.keys():
					pass	
				else:
					self.outer.__dict__[Name]=Function.__get__(self)
					self.outer.__methode_att__=[]	
					pass
			def Delete(Name="OMARI"):
					pass
		return Fun1(self)


def ImportFromFile(filename="19p.xlsx"):
	from xlrd import open_workbook
	wb = open_workbook(filename)
	sheet=wb.sheet_by_index(0)
	rowdata=[]
	for col in range(sheet.ncols):
		rowdata.append(sheet.col_values(colx=col))
		pass
	return rowdata

	
	

if __name__=="__main__":
	A=Attibute()
	A.value=[x for x in range(100)]
	print(A.Filter(lambda x:x>1).Return())
	print(A.Filter(lambda x:50>x>5).Return())
	print(A.value)
	A.Reset()
	print(A.value)
	B=ImportFromFile()
	A=database()
	def Fun(self):
		print ("OMARI -->")
	#A.Functions().Add(Fun,"Omari")
	