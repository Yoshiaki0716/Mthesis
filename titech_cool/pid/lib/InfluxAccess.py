from influxdb import InfluxDBClient
import pandas

class InfluxAccess:
    def __init__( self, config ):
        
        host        = config.get('host')
        port        = config.get('port')
        database    = config.get('database')
        measurement = config.get('measurement')
        if "tags" in config:
            self.tags = config.get('tags')
        else:
            self.tags = {}
        field       = config.get('field')
        
        self.client = InfluxDBClient( host = host, port = port, database = database )
        self.measurement = measurement
        self.field = field if field is not None else "*"

    def getLast( self ):
        try:
            if ' ' in self.field:
                field = '"{}"'.format( self.field )
            else:
                field = self.field
        except:
            pass
        if self.tags:
            where_parts = [f"{key}='{value}'" for key, value in self.tags.items()]
            where_clause = 'WHERE ' + ' AND '.join(where_parts)
        else:
            where_clause = ''
            query = 'SELECT {} FROM {} {} ORDER BY DESC LIMIT 1'.format( field, self.measurement, where_clause )
        result = self.client.query( query )
        df = pandas.DataFrame( result.get_points( measurement = self.measurement ) )
        self.last_data = df.values[-1]
        return float( df.values[-1][-1] )


    def getLastTime( self ):
        try:
            
            if ' ' in self.field:
                field = '"{}"'.format( self.field )
            else:
                field = self.field
        except:
            pass
        
        if self.tags:
            where_parts = [f'{key}="{value}"' for key, value in self.tags.items()]
            where_clause = 'WHERE ' + ' AND '.join(where_parts)
        else:
            where_clause = ''
        query = 'SELECT {} FROM {} {} ORDER BY DESC LIMIT 1'.format( field, self.measurement, where_clause )
        result = self.client.query( query )
        df = pandas.DataFrame( result.get_points( measurement = self.measurement ) )
        self.last_data = df.values[-1]
        return df.values[-1][0]


    def record( self, value ):
        
        data_ = [ {'fields': { self.field : value },
                   'measurement': self.measurement } ]

        res = self.client.write_points(data_)

