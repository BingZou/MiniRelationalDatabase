# Shiqing Li (sl7085)
# Bing Zou (bz1031)

import re
import os
import sys
import numpy as np
from BTrees.OOBTree import OOBTree
import copy
import time

class myDB:
    def __init__(self):
        self.inter_table = None 
        self.tablef_to_sort = None
        self.tables_to_sort = None
        self.arithop = ['+', '-', '*', '/']
        self.relop = ['=','!','>','<']
        self.flag = None
        self.index = {}
        self.hash = None
        self.btree = None
        self.vals =[]

    def inputfromfile(self, condition):
        '''
         (i)   This function find text file named from user specified input and create a np array table
         (ii)  Inputs: user should specify a name of the file and make sure it is in the current directory
         (iii) Outputs: output an array table with each record of a tuple
         (iv)  There is a global effect out of this function, a data array table is created 
        '''
        condition = condition.replace(' ','')
        condition = re.split('[,\(\)\n\s]',condition)
        unwanted = ['', 'inputfromfile','.txt','./'] # modified here 
        condition = [e for e in condition if e not in unwanted] 
        location = './'+condition[0]+'.txt'
#         location = os.path.expanduser(location)
        if os.path.exists(location):
            p = location
            with open(p, 'r') as f:
                firstLine = f.readline().strip('\n').split('|')
                col_num = len(firstLine)
                row_num = sum(1 for line in f)
                dt = [(i, '<i4') for i in firstLine[:(col_num-1)]]
                dt.append((firstLine[-1], 'U12'))
                arr = np.empty([0,0], dtype=dt)
            arr_list = []
            with open(p, 'r') as f:    
                for line in f.readlines()[1:]:
                    ls_line = line.strip('\n').split('|')
                    ls_line[:6] = [int(i) for i in ls_line[:6]]
                    ls_tp = tuple(i for i in ls_line)
                    arr_list.append(ls_tp)
            return np.asarray(arr_list, dtype = arr.dtype)
        else:
            print('invalid location')
        
    def select(self, condition):
        '''
         (i)   This function selects instances that matches the query input(s)
         (ii)  Inputs: inputs should specify which table instances will be selected from,
                and should have one or multiple queries
         (iii) Outputs: the output instances should be from the specified table from input and satisfy the conditions
                 if there are multiple queries with 'and', then instances return will match all conditions
                 or if 'or' in the queries, then instances that match any query will be returned
         (iv)  There is no side effect on the global scope. 
        '''

        flip = {'>=':'<=','<=':'>=','>':'<','<':'>'}
        org_condition = condition
        condition = condition.replace(' ','')
        condition = re.split('[,\(\)\n\s]',condition)        
        unwanted = ['', 'select']
        condition = [e for e in condition if e not in unwanted]
        if self.flag:
            if any(self.index[self.flag] in con for con in condition):
                self.__selectByHashOrIndex(org_condition)
                
        table = eval(condition[0])
        self.inter_table = table.copy()

        t_des = table.dtype.descr
        t_att = [x[0] for x in t_des]

        con_with_opr = [con for con in condition if any(x in con for x in self.arithop)]
        for con in con_with_opr:
            split_terms = re.split(r'(\W+)', con)
            attr =  list(set(split_terms).intersection(set(t_att)))[0]
            opr_idx = [i for i, opr in enumerate(split_terms) if opr in self.arithop][0]
            self.inter_table[attr] = eval('self.inter_table[attr]'+split_terms[opr_idx]+'int(split_terms[opr_idx+1])')  

        multi = False
        if 'or' in condition:
            mul_opr = '|'
            multi = True
        elif 'and' in condition: 
            mul_opr = '&'
            multi = True

        query = "self.inter_table["
        for i in range(len(condition)//2):
            j = 2*i+1
            sub_condition = re.split(r'(\W+)',condition[j])

            if any(x in sub_condition for x in self.arithop):
                ari = [i for i, sub in enumerate(sub_condition) if sub in self.arithop][0]
                sub_condition.remove(sub_condition[ari]) # remove arithop
                sub_condition.remove(sub_condition[ari]) # remove the constant after the arithop

            attr = list(set(sub_condition).intersection(set(t_att)))[0]
            opr = [sub for sub in sub_condition if any(x in sub for x in self.relop)][0]
            attr_idx = sub_condition.index(attr)
            opr_idx = sub_condition.index(opr)
            if attr_idx > opr_idx:
                opr = flip[opr]
                sub_condition.remove(flip[opr])
            else: sub_condition.remove(opr)
            sub_condition.remove(attr)

            if multi: query += '('

            if opr != '=':
                query += "self.inter_table['" + attr + "'" + "]" + opr + sub_condition[0] 
            else:
                query += "self.inter_table['" + attr + "'" + "]" + opr*2 + sub_condition[0]

            if multi: query += ')'+ mul_opr

        if multi:
            query = query[:-1] 
        query += "]"
        return eval(query)
                 
    def project(self, condition):
        ''' 
         (i)   This function selects all instances with attribute(s) from a table specified.
         (ii)  Inputs: inputs should specify which table instances will be selected from and interested attributes.
         (iii) Outputs: A table with all instances from the table, with selected attributes.
         (iv)  There no global effect on the other tables, but a new table is created. 
        '''
        condition = condition.replace(' ','')
        condition = re.split('[,\(\)\n\s]',condition)
        unwanted = ['', 'project']
        condition = [e for e in condition if e not in unwanted]
        table = eval(condition[0])        
        query = "table[["
        for i in range(len(condition)-1):
            j = i+1
            query = query +\
            "'" +\
            condition[j] +\
            "',"
        query = query[:-1] + "]]"    
        return eval(query)    
    
    
    def avg(self, condition):
        '''
         (i)   Given a table and attribute, this function computes the average of that attribute.
         (ii)  Inputs: inputs should specify which table and an attribute from the table user would like to compute average of.
         (iii) Outputs: output a scalar of the average of that attribute from a table
         (iv)  There is no global effect. 
        '''
        condition = condition.replace(' ','')
        condition = re.split('[,\(\)\n\s]',condition)
        unwanted = ['', 'avg']
        condition = [e for e in condition if e not in unwanted]
        table = eval(condition[0])
        avgRes = np.mean(table[condition[1]])
        return avgRes
        
    def sum_(self, condition):
        '''
         (i)   Given a table and attribute, this function computes the sum of that attribute.
         (ii)  Inputs: inputs should specify which table and attribute from the table user would like to compute sum of.
         (iii) Outputs: output a scalar of the sum of that attribute from a table
         (iv)  There is no global effect. 
        '''
        condition = condition.replace(' ','')
        condition = re.split('[,\(\)\n\s]',condition)
        unwanted = ['', 'sum']
        condition = [e for e in condition if e not in unwanted]
        table = eval(condition[0])
        sumRes = np.sum(table[condition[1]])
        return sumRes

    def count(self, condition):
        '''
         (i)   Given a table  this function counts the number of rows in that table
         (ii)  Inputs: inputs should specify which table the user would like to compute count from.
         (iii) Outputs: output a scalar which is the size of the table
         (iv)  There is no global effect. 
        '''
        condition = condition.replace(' ', '')
        condition = re.split('[,\(\)\n\s]', condition)
        unwanted = ['', 'count']
        condition = [e for e in condition if e not in unwanted]
        table = eval(condition[0])
        countRes = len(table[condition[1]])
        return countRes 
    
    def sumgroup(self, condition):
        '''
         (i)   Given a table and a list of attributes, this func computes the sum of the first attributes grouped by
         the rest of attributes.
         (ii)  Inputs: inputs should specify which table and attribute(s) from the table user would like to compute
             sum from, note that sum will be performed to the first attributes group by the rest of attributes
         (iii) Outputs: output a table that will have same number of attributes as the input attributes.
         (iv)  There is no global effect on other tables, but a table is created. 
        '''
        condition = condition.replace(' ','')
        condition = re.split('[,\(\)\n\s]',condition)
        unwanted = ['', 'sumgroup']
        condition = [e for e in condition if e not in unwanted]
        table = eval(condition[0])
        
        unique_list = []
        for con in condition[2:]:
            unique = np.unique(table[con])
            unique_list.append(unique)
            
        res_list = []
        
        if len(condition) == 3:
            unique = np.unique(table[condition[2]])
            for j in range(len(unique)):
                subset = table[table[condition[2]] == unique[j]]
                if subset.size > 0:
                    res_list.append(tuple((np.sum(subset[condition[1]]), unique[j])))
            dt = [('sum', 'i4'), (condition[2], table[condition[2]].dtype.str)]
                
        elif len(condition) == 4:
            dt = [('sum', 'i4'), (condition[2], table[condition[2]].dtype.str), (condition[3], table[condition[3]].dtype.str)]
            unique_first = np.unique(table[condition[2]])
            unique_second = np.unique(table[condition[3]])
            for i in unique_first:
                for j in unique_second:
                    subset = table[(table[condition[2]]==i)&(table[condition[3]]==j)]
                    if subset.size > 0:
                        res_list.append(tuple((np.sum(subset[condition[1]]), i, j)))
        return np.asarray(res_list, dtype = dt)
        
    def avggroup(self, condition):
        '''
         (i)   Given a table and a list of attributes, this func computes the average of the first attributes grouped by
         the rest of attributes.
         (ii)  Inputs: inputs should specify which table and attribute(s) from the table user would like to compute
             average from, note that average will be performed to the first attributes group by the rest of attributes
         (iii) Outputs: output a table that will have same number of attributes as the input attributes.
         (iv)  There is not global effect on other tables, but a new table is created. 
         '''
        
        condition = condition.replace(' ','')
        condition = re.split('[,\(\)\n\s]',condition)
        unwanted = ['', 'avggroup']
        condition = [e for e in condition if e not in unwanted]
        table = eval(condition[0])
        
        unique_list = []
        for con in condition[2:]:
            unique = np.unique(table[con])
            unique_list.append(unique)
            
        res_list = []
        
        if len(condition) == 3:
            unique = np.unique(table[condition[2]])
            for j in range(len(unique)):
                subset = table[table[condition[2]] == unique[j]]
                if subset.size > 0:
                    res_list.append(tuple((np.mean(subset[condition[1]]), unique[j])))
            dt = [('avg', 'f4'), (condition[2], table[condition[2]].dtype.str)]
                
        elif len(condition) == 4:
            dt = [('avg', 'f4'), (condition[2], table[condition[2]].dtype.str), (condition[3], table[condition[3]].dtype.str)]
            unique_first = np.unique(table[condition[2]])
            unique_second = np.unique(table[condition[3]])
            for i in unique_first:
                for j in unique_second:
                    subset = table[(table[condition[2]]==i)&(table[condition[3]]==j)]
                    if subset.size > 0:
                        res_list.append(tuple((np.mean(subset[condition[1]]), i, j)))
        return np.asarray(res_list, dtype = dt)
        
 
    
    def countgroup(self, condition):
        '''
         (i)   Given a table and a list of attributes, this func counts the number of first attributes grouped by
         the rest of attributes.
         (ii)  Inputs: inputs should specify which table and attribute(s) from the table user would like to compute
             count from, note that count will be performed to the first attributes group by the rest of attributes.
         (iii) Outputs: output a table that will have same number of attributes as the input attributes.
         (iv)  There is not global effect.
        '''
        condition = condition.replace(' ','')
        condition = re.split('[,\(\)\n\s]',condition)
        unwanted = ['', 'countgroup']
        condition = [e for e in condition if e not in unwanted]
        table = eval(condition[0])
        
        unique_list = []
        for con in condition[2:]:
            unique = np.unique(table[con])
            unique_list.append(unique)
            
        res_list = []
        
        if len(condition) == 3:
            unique = np.unique(table[condition[2]])
            for j in range(len(unique)):
                subset = table[table[condition[2]] == unique[j]]
                if subset.size > 0:
                    res_list.append(tuple((count, unique[j])))
            dt = [('count', 'i4'), (condition[2], table[condition[2]].dtype.str)]
                
        elif len(condition) == 4:
            dt = [('count', 'i4'), (condition[2], table[condition[2]].dtype.str), (condition[3], table[condition[3]].dtype.str)]
            unique_first = np.unique(table[condition[2]])
            unique_second = np.unique(table[condition[3]])
            for i in unique_first:
                for j in unique_second:
                    subset = table[(table[condition[2]]==i)&(table[condition[3]]==j)]
                    if subset.size > 0:
                        count = subset.size
                        res_list.append(tuple((count, i, j)))
        return np.asarray(res_list, dtype = dt)
    
    
    
    def join(self, condition):
        '''
         (i)   Given two tables and conditions(s), this func returns a joined table that matches all conditions. This
         function first sort both tables by the attributes that will perform condition with each time, and then check
         and return all instances where the condition is met. 
         (ii)  Inputs: inputs should specify which tables and condition(s) user would like to perform over.
         (iii) Outputs: output a joined table that each instance is a tuple with both tables information 
         (iv)  There is a global effect that a joined table is created. 
         '''
        if self.flag:
            self.joinByHashOrIndex(condition)
           
        condition = condition.replace(' ','')
        condition = re.split('[,\(\)\n\s]', condition)
        unwanted = ['', 'join']
        condition = [e for e in condition if e not in unwanted]
        table1 = eval(condition[0])
        table2 = eval(condition[1])
        self.tablef_to_sort = table1.copy()
        self.tables_to_sort = table2.copy()
        
        dt1 = self.tablef_to_sort.dtype.descr
        dt2 = self.tables_to_sort.dtype.descr
        dt1_sorted = [(condition[0] + '_' + x[0], x[1]) for x in dt1]
        dt2_sorted = [(condition[1] + '_'+ x[0], x[1]) for x in dt2]
        start_time = time.time()
        refined_condition = []
        for ind, con in enumerate(condition[2:]):
            if ind % 2 == 0: 
                sub_condition = re.split('[.\=\!=\<\>\<=\>=]', con)
                sub_condition = [sub for sub in sub_condition if sub != '']
                if any(x in sub_condition[1] for x in self.arithop):
                    ari = re.findall(r'[\+\-\*\/]', sub_condition[1])[0]
                    sub_condition[1], const = re.split(r'['+ari+']', sub_condition[1])
                    self.tablef_to_sort[sub_condition[1]] = eval('self.tablef_to_sort[sub_condition[1]]'+ ari +'int(const)')
                if any(x in sub_condition[3] for x in self.arithop):
                    ari = re.findall(r'[\+\-\*\/]', sub_condition[3])[0]
                    sub_condition[3], const = re.split(r'['+ari+']', sub_condition[3])
                    self.tables_to_sort[sub_condition[3]] = eval('self.tables_to_sort[sub_condition[3]]'+ ari +'int(const)')
                operation = re.findall(r'[\=\!\>\<]', con)
                operation = ''.join(operation)
                refined_condition.append([sub_condition[0], sub_condition[1], operation, sub_condition[2], sub_condition[3] ])
        
        self.tablef_to_sort = self.tablef_to_sort.astype(dt1_sorted)
        self.tables_to_sort = self.tables_to_sort.astype(dt2_sorted)

        dt_sorted = dt1_sorted + dt2_sorted
        query = 'self.tables_to_sort['
        table2_lookup_list = []
        table1_lookup_list = []
        for i, ref_con in enumerate(refined_condition):
            # table 1, col1, opr, table2, col2
            ref_con[1] = ref_con[0] +'_'+ref_con[1]
            ref_con[4] = ref_con[3] +'_'+ref_con[4]
            table2_lookup_list.append(ref_con[4])
            table1_lookup_list.append(ref_con[1])
        if ref_con[2] == '=':
            query += '(self.tables_to_sort["'+ref_con[4]+'"] == vals['+str(i)+'])&'
        else:
            query += '(self.tables_to_sort["'+ref_con[4]+'"] + ref_con[2]+ vals['+str(i)+'])&'
        query = query[:-1]
        query += ']'
        res_list = []
        for i in range(len(self.tablef_to_sort)):
            vals = [self.tablef_to_sort[j][i] for j in table1_lookup_list]
            table2_sel = eval(query)
            if table2_sel.size > 0:
                table2_sel = table2_sel.tolist()
                table1_repeat = [self.tablef_to_sort[i]] * len(table2_sel)
                toi = [tuple(table1_repeat[k])+tuple(table2_sel[k]) for k in range(len(table1_repeat))]
                res_list += toi
        return np.asarray(res_list, dtype = dt_sorted)
 
        
    def movavg(self, condition):
        '''
         (i)   Performs a moving average of an attribute from a table specified, the sliding window is user defined.
         (ii)  Inputs: Given an attribute from a table, and a sliding window size. 
         (iii) Outputs: output a table that performs moving average over an attribute. 
         (iv)  There is no global effect on other table, but a new table is created. 
        '''
        condition = condition.replace(' ','')
        condition = re.split('[,\(\)\n\s]',condition)
        unwanted = ['', 'movavg']
        condition = [e for e in condition if e not in unwanted]
        table = eval(condition[0])
        target_col = table[condition[1]]
        movAvg = []
        totalSum = 0
        count = eval(condition[2])
        for i in range(len(table)):
            if i < (count-1):
                totalSum += target_col[i]
                movAvg.append(totalSum/(i+1))
            else:
                totalSum = 0
                for j in range(count):
                    totalSum += target_col[i-j]
                movAvg.append(totalSum/count)
        
        dt = table.dtype.descr + [('move_avg', 'f4')]
        res_list = []
        for i in range(len(table)):
            toi = tuple(table[i])+(movAvg[i],)
            res_list.append(toi)
        return np.asarray(res_list, dtype = dt)
    
    def movsum(self, condition):
        '''
         (i)   Performs a moving sum of an attribute from a table specified, the sliding window is user defined.
         (ii)  Inputs: Given an attribute from a table, and a sliding window size. 
         (iii) Outputs: output a table that performs moving sum over an attribute. 
         (iv)  There is no global effect on other table, but a new table is created. 
        '''
        condition = condition.replace(' ','')
        condition = re.split('[,\(\)\n\s]',condition)
        unwanted = ['', 'movsum']
        condition = [e for e in condition if e not in unwanted]
        table = eval(condition[0])
        target_col = table[condition[1]]
        movSum = []
        totalSum = 0
        count = eval(condition[2])
        for i in range(len(table)):
            if i < (count-1):
                totalSum += target_col[i]
                movSum.append(totalSum)
            else:
                totalSum = 0
                for j in range(count):
                    totalSum += target_col[i-j]
                movSum.append(totalSum)
        dt = table.dtype.descr + [('move_sum', 'f4')]
        res_list = []
        for i in range(len(table)):
            toi = tuple(table[i])+(movSum[i],)
            res_list.append(toi)
        return np.asarray(res_list, dtype = dt)
    
    def Btree(self, condition, return_ = False):
        '''
         (i)   Create a btree with an index on an attribute from a table that user specified. 
         (ii)  Inputs: An attribute from a table
         (iii) Outputs: output a btree that is indexed with the attribute specified.
         (iv)  There is a global effect that a btree is created
        '''
        condition = condition.replace(' ','')
        condition = re.split('[,\(\)\n\s]',condition)
        unwanted = ['' ,'Btree']
        condition = [e for e in condition if e not in unwanted]
        table = eval(condition[0])
        btree = OOBTree()
        for i, value in enumerate(table[condition[1]]):
            if value in btree:
                btree.update({value: btree[value] + [i]})
            else:
                btree.update({value: [i]})
        self.flag = 'btree'
        self.index['btree'] = condition[1]
        self.btree = btree
        if return_:
            return btree
    
    def Hash(self, condition, return_ = False):
        '''
         (i)   Create a Hash with an index on an attribute from a table that user specified. 
         (ii)  Inputs: An attribute from a table
         (iii) Outputs: output a Hash that is indexed with the attribute specified.
         (iv)  There is a global effect that a hash table is created
        '''
        condition = condition.replace(' ','')
        condition = re.split('[,\(\)\n\s]',condition)
        unwanted = ['', 'hash','Hash']
        condition = [e for e in condition if e not in unwanted]
        table = eval(condition[0])        
        dic = {}
        for i, value in enumerate(table[condition[1]]):
            if value in dic:
                dic[value] += [i]
            else:
                dic[value] = [i]
        self.flag = 'hash'
        self.index['hash'] = condition[1]
        self.hash = dic
        if return_:
            return dic
        
    def __selectByHashOrIndex(self, condition):
        '''
         (i)   With already constructed Btree with index on a certain attribute(s), perform equality selection. 
         (ii)  Inputs: With a table constructed by Btree and an equality selection. 
         (iii) Outputs: output a table that matches the equality selection 
         (iv)  There is a global effect that a table is created. 
        '''
        
        condition = condition.replace(' ','')
        condition = re.split('[,\(\)\n\s]',condition)
        unwanted = ['', 'select']
        condition = [e for e in condition if e not in unwanted]
        sub_condition = re.split(r'(\W+)',condition[1])
        if sub_condition[1] == '=':
            #selectedIndex = self.btree[eval(sub_condition[2])]
            selectedIndex = eval('self.'+str(self.flag)+'['+sub_condition[2]+']')
            table = eval(condition[0])
            return table[selectedIndex]
        else:
            print("Invalid condition")
    

    def concat(self, condition):
        '''
         (i)   Given two tables return a joined table if the schemas from both tables are consistent.
         (ii)  Inputs: two tables that should have consistent schema.
         (iii) Outputs: output a concatenated table if two tables are consistent.
         (iv)  There is a global effect that a concat table is created. 
         '''
        condition = condition.replace(' ','')
        condition = re.split('[,\(\)\n\s]',condition)
        unwanted = ['', 'concat']
        condition = [e for e in condition if e not in unwanted]
        t1 = eval(condition[0])
        t2 = eval(condition[1]) 
        table1_schema = t1.dtype.descr
        table2_schema = t2.dtype.descr
        if any(s1 not in table2_schema for s1 in table1_schema) or any(s2 not in table1_schema for s2 in table2_schema):
            return print('inconsistent schema, unable to concat')
        return np.unique(np.concatenate((t1, t2), axis=0),0)
    
    
    def outputtofile(self, condition):
        '''
         (i)   Given a file path and a table, this function outputs the table into a file with "|" as seperator
         (ii)  Inputs: a file path and a table
         (iii) Outputs: nothing
         (iv)  There is a global effect a file is created locally 
         '''
        condition = condition.replace(' ','')
        condition = re.split('[,\(\)\n\s]',condition)
        unwanted = ['', 'outputtofile']
        condition = [e for e in condition if e not in unwanted]
        filename = condition[0] + '.txt'
        table=eval(condition[0])
        header = "|".join(str(x) for x in table.dtype.names)
        descr = table.dtype.descr
        fmt = []
        descr = table.dtype.descr
        for i in range(len(descr)):
            if descr[i][1] == '<i4': ### ?????
                fmt.append('%i')
            elif descr[i][1] == '<U12':
                fmt.append('%s')
        np.savetxt(filename, table,  fmt=fmt, delimiter='|', header=header, comments='')    
        
    def sort(self,condition):
        '''
        (i)   Sort returns a table that is sorted by attributes(s) users specified.
        (ii)  Inputs: inputs should specify which table instances will be sorted from,
                and should have one or multiple attributes and should be ordered by priority
        (iii) Outputs: the output schema should be from the specified table from input and 
                return a sorted table. 
        (iv)  There is no side effect on the global scope. 
        '''
        condition = condition.replace(' ','')
        condition = re.split('[,\(\)\n\s]',condition)
        unwanted = ['', 'sort']
        condition = [e for e in condition if e not in unwanted]

        table1 = eval(condition[0])
        table_sort = table1.copy()
        condition.remove(condition[0])

        return np.sort(table_sort, order = condition)  
    
    def joinByHashOrIndex(self, condition):
        '''
         (i)   Given two tables and a conditions, this func returns a joined table that matches the condition by 
                 Hash or Index, specify by previous operation
         (ii)  Inputs: inputs should specify which tables and condition user would like to perform over.
         (iii) Outputs: output a joined table that each instance is a tuple with both tables information 
         (iv)  There is a global effect that a joined table is created. 
         '''
        condition = condition.replace(' ','')
        condition = re.split('[,\(\)\n\s]',condition)
        unwanted = ['', 'join']
        condition = [e for e in condition if e not in unwanted]
        table1 = eval(condition[0])
        table2 = eval(condition[1])
        
        sub_condition = re.split('[.\=\!=\<\>\<=\>=]', condition[2])

        # dtypes
        dt1,dt2 = table1.dtype.descr,table2.dtype.descr
        dt1 = [(condition[0] + '_' + x[0], x[1]) for x in dt1]
        dt2 = [(condition[1] + '_'+ x[0], x[1]) for x in dt2]
        dt = dt1 + dt2
        
        res_list = []
        
        if self.flag == 'Btree': 
            for value in list(index_table1.keys()):
                if value in index_table2:
                    index_table1 = self.Btree(condition[0]+','+sub_condition[1], True)
                    index_table2 = self.Btree(condition[1]+','+sub_condition[3], True)         
                    indexes1 = index_table1[value]
                    indexes2 = index_table2[value]
                    tois = [tuple(table1[i1]) + tuple(table2[i2]) for i1 in indexes1 for i2 in indexes2]
                    res_list.append(tois)
                    
        elif self.flag == 'Hash':
            for value in index_table1.keys():
                if value in index_table2:
                    index_table1 = self.Hash(condition[0]+','+sub_condition[1], True)
                    index_table2 = self.Hash(condition[1]+','+sub_condition[3], True)  
                    indexes1 = index_table1[value]
                    indexes2 = index_table2[value]
                    tois = [tuple(table1[i1]) + tuple(table2[i2]) for i1 in indexes1 for i2 in indexes2]
                    res_list.append(tois)
                    
        return np.asarray(res_list, dtype = dt)
        

if __name__ == "__main__":
    # first initialize
    myclass=myDB()

    # read commands from file 
    input_commands = []
    with open("mycommand.txt", "r") as f:
        for line in f:
            line = line.strip().split("//")[0]
            if line:
                input_commands.append(line)

    # execute each command
    record_time, total_time_start = time.time(),time.time()
    for com in input_commands:
        if ':=' not in com:
            operation = re.split('[\(]',com.replace(' ',''))[0]
            exec('myclass.'+operation + '("' +com +'")', globals())
            print(f'execute {com} took {time.time()-record_time} secs')
            record_time = time.time()
        else:
            tablename, condition = com.split(":=")
            tablename = tablename.strip()
            operation = re.split('[,\(\)\n\s]',condition.replace(' ',''))
            operation = [op for op in operation if op != '']
            if operation[0] == 'sum':
                operation[0] = 'sum_'
            exec(tablename + "= myclass." + operation[0] + "(condition= '" + condition + "')", globals())
            print(f'execute {com} took {time.time()-record_time} secs')
            record_time = time.time()

    print(f'processed total {len(input_commands)} and took {time.time()-total_time_start} secs')

        