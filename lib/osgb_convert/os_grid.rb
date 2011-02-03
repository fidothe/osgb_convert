module OsgbConvert
  # /*
  #  * convert OS grid reference to geodesic co-ordinates
  #  */
  # function OSGridToLatLong(gridRef) {
  #   var gr = gridrefLetToNum(gridRef);
  #   var E = gr[0], N = gr[1];
  # 
  #   var a = 6377563.396, b = 6356256.910;              // Airy 1830 major & minor semi-axes
  #   var F0 = 0.9996012717;                             // NatGrid scale factor on central meridian
  #   var lat0 = 49*Math.PI/180, lon0 = -2*Math.PI/180;  // NatGrid true origin
  #   var N0 = -100000, E0 = 400000;                     // northing & easting of true origin, metres
  #   var e2 = 1 - (b*b)/(a*a);                          // eccentricity squared
  #   var n = (a-b)/(a+b), n2 = n*n, n3 = n*n*n;
  # 
  #   var lat=lat0, M=0;
  #   do {
  #     lat = (N-N0-M)/(a*F0) + lat;
  # 
  #     var Ma = (1 + n + (5/4)*n2 + (5/4)*n3) * (lat-lat0);
  #     var Mb = (3*n + 3*n*n + (21/8)*n3) * Math.sin(lat-lat0) * Math.cos(lat+lat0);
  #     var Mc = ((15/8)*n2 + (15/8)*n3) * Math.sin(2*(lat-lat0)) * Math.cos(2*(lat+lat0));
  #     var Md = (35/24)*n3 * Math.sin(3*(lat-lat0)) * Math.cos(3*(lat+lat0));
  #     M = b * F0 * (Ma - Mb + Mc - Md);                // meridional arc
  # 
  #   } while (N-N0-M >= 0.00001);  // ie until < 0.01mm
  # 
  #   var cosLat = Math.cos(lat), sinLat = Math.sin(lat);
  #   var nu = a*F0/Math.sqrt(1-e2*sinLat*sinLat);              // transverse radius of curvature
  #   var rho = a*F0*(1-e2)/Math.pow(1-e2*sinLat*sinLat, 1.5);  // meridional radius of curvature
  #   var eta2 = nu/rho-1;
  # 
  #   var tanLat = Math.tan(lat);
  #   var tan2lat = tanLat*tanLat, tan4lat = tan2lat*tan2lat, tan6lat = tan4lat*tan2lat;
  #   var secLat = 1/cosLat;
  #   var nu3 = nu*nu*nu, nu5 = nu3*nu*nu, nu7 = nu5*nu*nu;
  #   var VII = tanLat/(2*rho*nu);
  #   var VIII = tanLat/(24*rho*nu3)*(5+3*tan2lat+eta2-9*tan2lat*eta2);
  #   var IX = tanLat/(720*rho*nu5)*(61+90*tan2lat+45*tan4lat);
  #   var X = secLat/nu;
  #   var XI = secLat/(6*nu3)*(nu/rho+2*tan2lat);
  #   var XII = secLat/(120*nu5)*(5+28*tan2lat+24*tan4lat);
  #   var XIIA = secLat/(5040*nu7)*(61+662*tan2lat+1320*tan4lat+720*tan6lat);
  # 
  #   var dE = (E-E0), dE2 = dE*dE, dE3 = dE2*dE, dE4 = dE2*dE2, dE5 = dE3*dE2, dE6 = dE4*dE2, dE7 = dE5*dE2;
  #   lat = lat - VII*dE2 + VIII*dE4 - IX*dE6;
  #   var lon = lon0 + X*dE - XI*dE3 + XII*dE5 - XIIA*dE7;
  # 
  #   return new LatLon(lat.toDeg(), lon.toDeg());
  # }
  # 
  # 
  # /* 
  #  * convert standard grid reference ('SU387148') to fully numeric ref ([438700,114800])
  #  *   returned co-ordinates are in metres, centred on grid square for conversion to lat/long
  #  *
  #  *   note that northern-most grid squares will give 7-digit northings
  #  *   no error-checking is done on gridref (bad input will give bad results or NaN)
  #  */
  # function gridrefLetToNum(gridref) {
  #   // get numeric values of letter references, mapping A->0, B->1, C->2, etc:
  #   var l1 = gridref.toUpperCase().charCodeAt(0) - 'A'.charCodeAt(0);
  #   var l2 = gridref.toUpperCase().charCodeAt(1) - 'A'.charCodeAt(0);
  #   // shuffle down letters after 'I' since 'I' is not used in grid:
  #   if (l1 > 7) l1--;
  #   if (l2 > 7) l2--;
  # 
  #   // convert grid letters into 100km-square indexes from false origin (grid square SV):
  #   var e = ((l1-2)%5)*5 + (l2%5);
  #   var n = (19-Math.floor(l1/5)*5) - Math.floor(l2/5);
  # 
  #   // skip grid letters to get numeric part of ref, stripping any spaces:
  #   gridref = gridref.slice(2).replace(/ /g,'');
  # 
  #   // append numeric part of references to grid index:
  #   e += gridref.slice(0, gridref.length/2);
  #   n += gridref.slice(gridref.length/2);
  # 
  #   // normalise to 1m grid, rounding up to centre of grid square:
  #   switch (gridref.length) {
  #     case 6: e += '50'; n += '50'; break;
  #     case 8: e += '5'; n += '5'; break;
  #     // 10-digit refs are already 1m
  #   }
  # 
  #   return [e, n];
  # }
  # 
  # 
  # /*
  #  * convert numeric grid reference (in metres) to standard-form grid ref
  #  */
  # function gridrefNumToLet(e, n, digits) {
  #   // get the 100km-grid indices
  #   var e100k = Math.floor(e/100000), n100k = Math.floor(n/100000);
  # 
  #   if (e100k<0 || e100k>6 || n100k<0 || n100k>12) return '';
  # 
  #   // translate those into numeric equivalents of the grid letters
  #   var l1 = (19-n100k) - (19-n100k)%5 + Math.floor((e100k+10)/5);
  #   var l2 = (19-n100k)*5%25 + e100k%5;
  # 
  #   // compensate for skipped 'I' and build grid letter-pairs
  #   if (l1 > 7) l1++;
  #   if (l2 > 7) l2++;
  #   var letPair = String.fromCharCode(l1+'A'.charCodeAt(0), l2+'A'.charCodeAt(0));
  # 
  #   // strip 100km-grid indices from easting & northing, and reduce precision
  #   e = Math.floor((e%100000)/Math.pow(10,5-digits/2));
  #   n = Math.floor((n%100000)/Math.pow(10,5-digits/2));
  # 
  #   var gridRef = letPair + e.padLZ(digits/2) + n.padLZ(digits/2);
  # 
  #   return gridRef;
  # }
  class OSGrid
    GRID_LETTERS = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'J', 'K',
                    'L', 'M', 'N', 'O', 'P', 'Q', 'S', 'T', 'U', 'V', 
                    'W', 'X', 'Y', 'Z']

    attr_reader :easting, :northing

    def self.from_wgs84(wgs84)
      from_osgb36(wgs84.osgb36)
    end

    # http://www.movable-type.co.uk/scripts/latlong-gridref.html
    # (c) Chris Veness 2005-2010  Released under an LGPL license
    # http://www.fsf.org/licensing/licenses/lgpl.html
    # Ported to ruby by Harry Wood

    # OSGB36 coordinates to OS UK grid eastings & northings
    def self.from_osgb36(osgb36)
      lat = OsgbConvert.degrees_to_rads(osgb36.lat);
      lon = OsgbConvert.degrees_to_rads(osgb36.long);

      a = 6377563.396; b = 6356256.910          # Airy 1830 major & minor semi-axes
      f0 = 0.9996012717                         # NatGrid scale factor on central meridian
      lat0 = OsgbConvert.degrees_to_rads(49); lon0 = OsgbConvert.degrees_to_rads(-2)        # NatGrid true origin
      n0 = -100000; e0 = 400000;                # northing & easting of true origin, metres
      e2 = 1 - (b*b) / (a*a);                   # eccentricity squared
      n = (a-b) / (a+b); n2 = n*n; n3 = n*n*n;

      cos_lat = Math.cos(lat); sinLat = Math.sin(lat);
      nu = a*f0/Math.sqrt(1-e2*sinLat*sinLat);              # transverse radius of curvature
      rho = a*f0*(1-e2) / ( (1-e2*sinLat*sinLat) ** 1.5);  # meridional radius of curvature
      eta2 = nu/rho-1;

      ma = (1 + n + (5/4)*n2 + (5/4)*n3) * (lat-lat0);
      mb = (3*n + 3*n*n + (21/8)*n3) * Math.sin(lat-lat0) * Math.cos(lat+lat0);
      mc = ((15/8)*n2 + (15/8)*n3) * Math.sin(2*(lat-lat0)) * Math.cos(2*(lat+lat0));
      md = (35/24)*n3 * Math.sin(3*(lat-lat0)) * Math.cos(3*(lat+lat0));
      m = b * f0 * (ma - mb + mc - md);              # meridional arc

      cos3lat = cos_lat*cos_lat*cos_lat;
      cos5lat = cos3lat*cos_lat*cos_lat;
      tan2lat = Math.tan(lat)*Math.tan(lat);
      tan4lat = tan2lat*tan2lat;

      i = m + n0;
      ii = (nu/2)*sinLat*cos_lat;
      iii = (nu/24)*sinLat*cos3lat*(5-tan2lat+9*eta2);
      iiiA = (nu/720)*sinLat*cos5lat*(61-58*tan2lat+tan4lat);
      iv = nu*cos_lat;
      v = (nu/6)*cos3lat*(nu/rho-tan2lat);
      vi = (nu/120) * cos5lat * (5 - 18*tan2lat + tan4lat + 14*eta2 - 58*tan2lat*eta2);

      dLon = lon-lon0;
      dLon2 = dLon*dLon
      dLon3 = dLon2*dLon
      dLon4 = dLon3*dLon
      dLon5 = dLon4*dLon
      dLon6 = dLon5*dLon

      n = i + ii*dLon2 + iii*dLon4 + iiiA*dLon6;
      e = e0 + iv*dLon + v*dLon3 + vi*dLon5;

      new(e, n) # return new OSGrid instance using the raw easting and northings
    end

    # convert standard grid reference ('SU387148') to fully numeric ref ([438700,114800])
    # returned co-ordinates are in metres, centred on grid square for conversion to lat/long
    # 
    # note that northern-most grid squares will give 7-digit northings
    # no error-checking is done on gridref (bad input will give bad results or NaN)
    def self.from_standard_ref(grid_ref)
      # get numeric values of letter references, mapping A->0, B->1, C->2, etc:
      first_letter  = GRID_LETTERS.index(grid_ref[0..0].upcase)
      second_letter = GRID_LETTERS.index(grid_ref[1..1].upcase)

      # convert grid letters into 100km-square indexes from false origin (grid square SV):
      easting  = ((first_letter - 2) % 5) * 5 + (second_letter % 5)
      northing = (19 - (first_letter / 5).floor * 5) - (second_letter / 5).floor

      # skip grid letters to get numeric part of ref, stripping any spaces:
      numeric_ref = grid_ref[2..-1].gsub(/\s/, '')

      # append numeric part of references to grid index:
      easting  += numeric_ref[0..(numeric_ref.length / 2))
      northing += numeric_ref[(numeric_ref.length / 2)..-1]

      # normalise to 1m grid, rounding up to centre of grid square:
      case numeric_ref.length
      when 6
        easting  += '50'
        northing += '50'
      when 8
        easting  += '5'
        northing += '5'
      end
      # 10-digit refs are already 1m

      new(e, n)
    end

    def initialize(easting, northing)
      @easting = easting
      @northing = northing
    end

    # /*
    #  * convert OS grid reference to geodesic co-ordinates
    #  */
    # function OSGridToLatLong(gridRef) {
    # 
    #   tanLat = Math.tan(lat);
    #   tan2lat = tanLat*tanLat, tan4lat = tan2lat*tan2lat, tan6lat = tan4lat*tan2lat;
    #   secLat = 1/cosLat;
    #   nu3 = nu*nu*nu, nu5 = nu3*nu*nu, nu7 = nu5*nu*nu;
    #   VII = tanLat/(2*rho*nu);
    #   VIII = tanLat/(24*rho*nu3)*(5+3*tan2lat+eta2-9*tan2lat*eta2);
    #   IX = tanLat/(720*rho*nu5)*(61+90*tan2lat+45*tan4lat);
    #   X = secLat/nu;
    #   XI = secLat/(6*nu3)*(nu/rho+2*tan2lat);
    #   XII = secLat/(120*nu5)*(5+28*tan2lat+24*tan4lat);
    #   XIIA = secLat/(5040*nu7)*(61+662*tan2lat+1320*tan4lat+720*tan6lat);
    #   
    #   dE = (E-E0), dE2 = dE*dE, dE3 = dE2*dE, dE4 = dE2*dE2, dE5 = dE3*dE2, dE6 = dE4*dE2, dE7 = dE5*dE2;
    #   lat = lat - VII*dE2 + VIII*dE4 - IX*dE6;
    #   lon = lon0 + X*dE - XI*dE3 + XII*dE5 - XIIA*dE7;
    # 
    #   return new LatLon(lat.toDeg(), lon.toDeg());
    # }
    def osgb36
      # Airy 1830 major & minor semi-axes
      a = Converter::ELLIPSE[:airy1830][:a]
      b = Converter::ELLIPSE[:airy1830][:b]
      f0 = 0.9996012717;                      # NatGrid scale factor on central meridian
      lat0 = 49 * Math::PI / 180              # NatGrid true origin
      lon0 = -2 * Math::PI / 180
      n0 = -100000                            # northing of true origin, metres
      e0 = 400000                             # easting of true origin, metres
      e2 = 1 - (b * b) / (a * a)              # eccentricity squared
      n = (a - b) / (a + b)
      n2 = n * n, 
      n3 = n * n * n

      lat = lat0
      m = 0

      while (northing - n0 - m) >= 0.00001 # ie until < 0.01mm
        lat = (northing - n0 - m) / (a * f0) + lat

        var ma = (1 + n + (5 / 4) * n2 + (5 / 4) * n3) * (lat - lat0)
        var mb = (3 * n + 3 * n * n + (21 / 8) * n3) * Math.sin(lat - lat0) * Math.cos(lat + lat0)
        var mc = ((15 / 8) * n2 + (15 / 8) * n3) * Math.sin(2 * (lat - lat0)) * Math.cos(2 * (lat + lat0))
        var md = (35 / 24) * n3 * Math.sin(3 * (lat - lat0)) * Math.cos(3 * (lat + lat0))
        m = b * f0 * (ma - mb + mc - md) # meridional arc
      end

      cos_lat = Math.cos(lat), sin_lat = Math.sin(lat)
      nu = a * f0 / Math.sqrt(1 - e2 * sin_lat * sin_lat)           # transverse radius of curvature
      rho = a * f0 * (1 - e2) / (1 - e2 * sin_lat * sin_lat) ** 1.5 # meridional radius of curvature
      eta2 = nu / rho - 1

      tan_lat = Math.tan(lat);
      tan2lat = tan_lat*tan_lat
      tan4lat = tan2lat*tan2lat
      tan6lat = tan4lat*tan2lat
      sec_lat = 1/cos_lat
      nu3 = nu*nu*nu
      nu5 = nu3*nu*nu
      nu7 = nu5*nu*nu
      vii = tan_lat/(2*rho*nu)
      viii = tan_lat/(24*rho*nu3)*(5+3*tan2lat+eta2-9*tan2lat*eta2)
      ix = tan_lat/(720*rho*nu5)*(61+90*tan2lat+45*tan4lat)
      x = sec_lat/nu
      xi = sec_lat/(6*nu3)*(nu/rho+2*tan2lat)
      xii = sec_lat/(120*nu5)*(5+28*tan2lat+24*tan4lat)
      xiia = sec_lat/(5040*nu7)*(61+662*tan2lat+1320*tan4lat+720*tan6lat)

      dE = (easting - e0)
      dE2 = dE*dE
      dE3 = dE2*dE
      dE4 = dE2*dE2
      dE5 = dE3*dE2
      dE6 = dE4*dE2
      dE7 = dE5*dE2
      lat = OsgbConvert.rads_to_degrees(lat - vii*dE2 + viii*dE4 - ix*dE6)
      lon = OsgbConvert.rads_to_degrees(lon0 + x*dE - xi*dE3 + xii*dE5 - xiia*dE7)
      
      OSGB36.new(lat, lon, 0)
    end

    # convert numeric grid reference (in metres) to standard-form grid ref
    # [defaults to 8 digits (because this was in the example, not because I know better)]
    def grid_ref(digits = 8)
      return @grid_ref if @grid_ref
      #get the 100km-grid indices
      e100k = (easting / 100000).floor
      n100k = (northing / 100000).floor

      return '' if (e100k < 0 or e100k > 6 or n100k < 0 or n100k > 12)

      #translate those into numeric equivalents of the grid letters
      first_letter  = (19 - n100k) - (19 - n100k) % 5 + ((e100k + 10) / 5).floor;
      second_letter = (19 - n100k) * 5 % 25 + e100k % 5;

      # letter - 1 to ensure we have 0-indexed the array and aren't off-by-one
      grid_name =  GRID_LETTERS[first_letter - 1] + GRID_LETTERS[second_letter - 1]

      # strip 100km-grid indices from easting & northing, and reduce precision
      reduced_precision_easting = ( (easting  % 100000) / (10 ** (5 - digits / 2)) ).floor
      reduced_precision_northing = ( (northing % 100000) / (10 ** (5 - digits / 2)) ).floor

      @grid_ref = grid_name + reduced_precision_easting.to_s.rjust(digits / 2) + reduced_precision_northing.to_s.rjust(digits / 2)
    end
  end
end