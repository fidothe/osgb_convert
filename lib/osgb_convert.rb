module OsgbConvert
  def self.degrees_to_rads(degrees)
    degrees * Math::PI / 180
  end

  def self.rads_to_degrees(rads)
    rads * 180 / Math::PI 
  end

  module Converter
    #ellipse parameters
    ELLIPSE = { 
      :wgs84 =>    { :a=> 6378137,     :b=> 6356752.3142, :f=> 1 / 298.257223563 },
      :airy1830 => { :a=> 6377563.396, :b=> 6356256.910,  :f=> 1 / 299.3249646   } 
    }

    #helmert transform parameters
    HELMERT = { 
      :wgs84toOSGB36 => { 
        :tx=> -446.448,  :ty=>  125.157,  :tz=> -542.060,  # m
        :rx=>  -0.1502,  :ry=>   -0.2470, :rz=>  -0.8421,  # sec
        :s=>   20.4894                                     # ppm
      },
      :osgb36toWGS84 => { 
        :tx=>  446.448,  :ty=> -125.157,  :tz=>  542.060,
        :rx=>    0.1502, :ry=>    0.2470, :rz=>    0.8421,
        :s=>   -20.4894
      }
    }

    def convert(p1lat, p1lon, p1height, e1, t, e2)
       # -- convert polar to cartesian coordinates (using ellipse 1)

       p1lat = OsgbConvert.degrees_to_rads(p1lat); p1lon = OsgbConvert.degrees_to_rads(p1lon); 

       a = e1[:a]; b = e1[:b];

       sinPhi = Math.sin(p1lat); cosPhi = Math.cos(p1lat);
       sinLambda = Math.sin(p1lon); cosLambda = Math.cos(p1lon);
       h = p1height;

       eSq = (a*a - b*b) / (a*a);
       nu = a / Math.sqrt(1 - eSq*sinPhi*sinPhi);

       x1 = (nu+h) * cosPhi * cosLambda;
       y1 = (nu+h) * cosPhi * sinLambda;
       z1 = ((1-eSq)*nu + h) * sinPhi;

       # -- apply helmert transform using appropriate params

       tx = t[:tx]; ty = t[:ty]; tz = t[:tz];
       rx = t[:rx] / 3600 * Math::PI/180;  #normalise seconds to radians
       ry = t[:ry] / 3600 * Math::PI/180;
       rz = t[:rz] / 3600 * Math::PI/180;
       s1 = t[:s] / 1e6 + 1;              #normalise ppm to (s+1)

       #apply transform
       x2 = tx + x1*s1 - y1*rz + z1*ry;
       y2 = ty + x1*rz + y1*s1 - z1*rx;
       z2 = tz - x1*ry + y1*rx + z1*s1;

       # -- convert cartesian to polar coordinates (using ellipse 2)

       a = e2[:a]; b = e2[:b];
       precision = 4 / a;  # results accurate to around 4 metres

       eSq = (a*a - b*b) / (a*a);
       p = Math.sqrt(x2*x2 + y2*y2);
       phi = Math.atan2(z2, p*(1-eSq)); phiP = 2 * Math::PI;
       while ( (phi-phiP).abs > precision) do
          nu = a / Math.sqrt(1 - eSq*Math.sin(phi)*Math.sin(phi));
          phiP = phi;
          phi = Math.atan2(z2 + eSq*nu*Math.sin(phi), p);
       end
       lambda = Math.atan2(y2, x2);
       h = p/Math.cos(phi) - nu;

       #return array [lat,lon,height]
       return [ OsgbConvert.rads_to_degrees(phi), OsgbConvert.rads_to_degrees(lambda), h ]; 
    end
  end
  
  class Coordinate
    include Converter
    
    attr_reader :lat, :long, :height
    
    def initialize(lat, long, height)
      @lat = lat
      @long = long
      @height = height
    end
  end
  
  class WGS84 < Coordinate
    def osgb36
      OSGB36.new(*convert(lat, long, height, ELLIPSE[:wgs84], HELMERT[:wgs84toOSGB36], ELLIPSE[:airy1830]))
    end
  end
  
  class OSGB36 < Coordinate
    def wgs84
      WGS84.new(*convert(lat, long, height, ELLIPSE[:airy1830], HELMERT[:osgb36toWGS84], ELLIPSE[:wgs84]))
    end
  end
end

require 'osgb_convert/os_grid'