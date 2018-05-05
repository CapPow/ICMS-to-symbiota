
# coding: utf-8

# #### Note:
# 
# ICMS uses a field, "Other Numbers" which would be appropriate to list Symbiota catalog Number, # (ie:	"GSMNP00030")
# 
# Also vise versa is appropriate for symbiota "otherCatalogNumber" housing the ICMS catalog number (ie: "GRSM  102763")

# In[805]:


import pandas as pd
import numpy as np
import re  # numeric extraction
from datetime import datetime # date alignment
import requests # make request to Taxonomic Resolution service "http://tnrs.iplantcollaborative.org/api.html"
import json # decode that request


# In[806]:


df = pd.read_csv('ICMS Data Export Example 042518.csv', dtype = 'str')
df['Cat #'] = df['Catalog #'].str.replace('GRSM','') # strip out their 'GRSM'
df['Cat #'] = df['Cat #'].str.strip() # remove the variable quantity of whitespace
refList = pd.read_csv('GSMNP DBLink 4-27-2018.csv', dtype = 'str')
refList = refList.rename({'GSMNP (number only)':'Cat #'}, axis = 'columns')
dfM = pd.merge(refList, df, how='inner') # create a "Merged data frame"


# In[807]:


#colMapper is structured as: {ICMS:Symbiota}
colMapper = {
    'Catalog #':'otherCatalogNumber',
    'SERNEC':'catalogNumber',
    'Cataloger':'recordEnteredBy',
    'Catalog Date':'dateEntered', #needs date converted
    'Ident Date':'dateIdentified', #needs date converted
    'Collection Date':'eventDate', #needs date converted
    'Elevation':'minimumElevationInMeters', # needs cleaned and converted
    'Collector':'recordedBy', #needs cleaned WRT "et al" & Possibly lists of names
    'Assoc Spec':'associatedTaxa', #some of these are probably also mixed in with "Description"
    'County':'county',
    'State':'stateProvince', #needs converted from TN, NC .. etc .. etc..
    'Locality':'locality',
    'Habitat':'habitat',
    'Identified By':'identifiedBy',
    'Location':'disposition', # possibly subjective assignment, get second opinion.
    'Kingdom':'kingdom',
    'Sci. Name:Genus':'genus',
    'Sci. Name:Species':'specificEpithet',
    'Sci. Name':'scientificName',
    'Sci. Name:Species Authority':'scientificNameAuthorship',
    'Family':'family',
    'Sex':'sex'}

dfS = dfM.rename(colMapper, axis = 'columns') # create a Symbiota dataFrame, and rename the fields accordingly.


# Field cleaning functions

# In[808]:


# clean fields, require a cell value
def expandStateName(stateStr):
    """expects a state abbreviation, returns an entire state name"""
    stateStr = stateStr.upper()
    refDict = {'TN':'Tennessee','NC':'North Carolina'}
    return refDict.get(stateStr)

def elevationToMeters(elevationStr):
    """expects elevation in feet, accepts additional non-numerical characters,
    prepares the value for DWC: minimumElevationInMeters """
    elevationStr = str(elevationStr).lower() # make sure the elevation is a string
    number = re.findall(r'\d+', elevationStr) # extract the numbers from it
    if len(number) == 0:
        return np.nan # if no numbers were found, send a nan down the line.
    number = min([ int(x) if x.isdigit() else x for x in number ]) # find the lowest continuous number
    if 'm' not in elevationStr: # if there is no letter "m" assume it is in feet.
        eleInMeters = (0.3048 * number) # and convert feet to meters
    else: # if "m" is in it
        eleInMeters = number
        
    eleInMeters = round(eleInMeters,1)
    return eleInMeters

def dateConverter(dateObj):
    """accepts a date or 4 digit year and tries to return it as ISO YYYY-MM-DD, otherwise returns the unchanged value"""
    cleanedDate = str(dateObj).strip()
    if cleanedDate == 'nan':
        return dateObj
    if (len(cleanedDate) == 4) & (cleanedDate.isdigit()):
        return '{}-00-00'.format(cleanedDate)
    else:
        dateObj = pd.to_datetime(cleanedDate, errors = 'ignore')
        try:
            result = dateObj.strftime('%Y-%m-%d')
        except ValueError:
            result = cleanedDate
        return result
    
# generate combination fields require an entire row (series) of data
def gen_lifeStage(row):
    """expects row data delivered as series returns a generated value"""
    potCols = ['Age/Stage','sex'] # group of potentitally useful column names
    potCols = [row[x] for x in potCols] # convert column names to values at those columns    
    colList = [x for x in potCols if not pd.isna(x)] # only keep stuff which is not empty
    age = str(row['Age'])
    if (age != np.nan) & (age not in str(row['Age/Stage'])): # if age is also not empty AND not already represented in 'Age/Stage'
        colList.insert(0, 'age: {}. '.format(age)) # then add it to the front of the list.

    result = ', '.join(colList).strip() # join the list into a comma seperated string.
    return result

def gen_locality(row):
    """expects row data delivered as series returns a generated value"""
    potCols = ['I-M Network','stateProvince','county'] # group of potentitally useful column names
    potCols = [row[x] for x in potCols] # convert column names to values at those columns    
    colList = [x for x in potCols if not pd.isna(x)] # only keep stuff which is not empty
    
    locality = str(row['locality'])
    if (locality != np.nan) & (locality not in potCols): # if age is also not empty AND not already represented in 'Age/Stage'
        colList.append(locality) # then add it to the front of the list.
        
    result = ', '.join(colList).strip() # join the list into a comma seperated string.
    return result

def gen_sciName(row):
    """expects row data delivered as series returns a tuble of generated values"""
    
    potCols = ['genus','specificEpithet','Sci. Name:Subspecies','Sci. Name:Variety'] # group of potentitally useful column names
    potCols = [row[x] for x in potCols] # convert column names to values at those columns
    colList = [x for x in potCols if not pd.isna(x)] # only keep stuff which is not empty
    queryName = '%20'.join(colList)
    url = 'http://tnrs.iplantc.org/tnrsm-svc/matchNames?retrieve=best&names={}'.format(queryName)
    result = json.loads(requests.get(url=url).content).get('items')[0]
    sciName = result.get('acceptedName')
    sciAuthor = result.get('acceptedAuthor')
    score = float(result.get('scientificScore')) # the confidence in the return
    if score < .98: # warn about issues with the names
        print('Scientific name for ICMS Catalog #: {}, Symbiota Catalog #: {}, scored: {}'.format(row['otherCatalogNumber'],row['catalogNumber'],score))
    result = (sciName, sciAuthor)
    return result

def gen_GPS(row):
    """expects row delivered as series returns a group of generated values"""
    
    if not pd.isna(row['UTM Z/E/N']): # if they used UTM coords, place them into verbatum coords
        verbCoord = 'UTM: {}'.format(row['UTM Z/E/N'])
        latLon = [np.nan,np.nan]
    if not pd.isna(row['Lat LongN/W']): # Check all the various coord fields for which ones are used
        verbCoord = row['Lat LongN/W']
        latLon = [float(x.upper().replace('N','').replace('W','')) for x in verbCoord.split('/')]
    elif not pd.isna(row['Lat LongN/W:Latitude Longitude']):
        verbCoord = row['Lat LongN/W:Latitude Longitude']
        latLon = [float(x.upper().replace('N','').replace('W','')) for x in verbCoord.split('/')]
    elif not pd.isna([row['Lat LongN/W:Latitude Degree'],row['Lat LongN/W:Longitude Degree']]).all(): # if deg / min/ sec used, handle it.
        latList = ['Lat LongN/W:Latitude Degree','Lat LongN/W:Latitude Minutes','Lat LongN/W:Latitude Seconds']
                                 #"[^0-9^.]"
        a, b, c = [float(re.sub('[^0-9^.]', '',row[i])) if re.sub('[^0-9^.]', '',str(row[i])).isnumeric() else 0 for i in latList]
        lat = a + (b / 60) + (c / 3600)     # make the conversions
        av,bv,cv = [row[i] if not pd.isna(row[i]) else '' for i in latList] # Keep the unmodified verbatim values
        lonList = ['Lat LongN/W:Longitude Degree','Lat LongN/W:Longitude Minutes','Lat LongN/W:Longitude Seconds']
        x, y, z = [float(re.sub('[^0-9^.]', '',row[i])) if re.sub('[^0-9^.]', '',str(row[i])).isnumeric() else 0 for i in lonList]
        xv,yv,zv = [row[i] if not pd.isna(row[i]) else '' for i in lonList] # Keep the unmodified verbatim values
        lon = x + (y / 60) + (z / 3600)     # make the conversions 
        verbCoord = 'lat: {}°, {}’, {}” lon:{}°, {}’, {}”'.format(av,bv,cv,xv,yv,zv) # format verbatim coords
        latLon = [lat,lon]
    else:
        verbCoord = np.nan # if all else fails, return nans
        latLon = [np.nan,np.nan]

    lat, lon = latLon
    if lon > 0: # since these are all national parks, safe to say it'll be in the western Hemisphere
        lon = -lon
    result = (lat,lon,verbCoord)
    return result


# In[809]:


# clean up all the date fields
for dateCol in ['dateEntered','dateIdentified','eventDate']:
    dfS[dateCol] = dfS[dateCol].apply(dateConverter)

# clean up the elevation
dfS['minimumElevationInMeters'] = dfS['minimumElevationInMeters'].apply(elevationToMeters)

# clean up the state column (expand the abbreviation)
dfS['stateProvince'] = dfS['stateProvince'].apply(expandStateName)

# clean up the county field
dfS['county'] = dfS['county'].str.title()

# generate the lifeStage column
dfS['lifeStage'] = dfS.apply(gen_lifeStage, axis=1)

# generate the Locality column
dfS['locality'] = dfS.apply(gen_locality, axis=1)

# generate the scientificName column
dfS[['scientificName','scientificNameAuthorship']] = dfS.apply(lambda row: pd.Series(gen_sciName(row)), axis=1)

# generate the gps coords & preserve the verbatim as best as possible.
dfS[['decimalLatitude', 'decimalLongitude', 'verbatimCoordinates']] = dfS.apply(lambda row: pd.Series(gen_GPS(row)), axis=1)


# In[810]:


colList = list(colMapper.values()) + ['lifeStage','decimalLatitude', 'decimalLongitude', 'verbatimCoordinates'] # select desired columns
dfS = dfS.loc[:,colList].copy() # reduce the Dataframe to that with Symbiota fields
print(dfS)


# In[812]:


dfS.to_csv('output.csv', encoding = 'utf-8')

