
import Foundation 
import MathTools

// This module supports fast neighbor lookup for atoms or other objects characterized by location and radius

// Essential ideas :
//
//      1) each atom/object has coordinate and radius ; find min/max for each coordinate
//      2) have a spatial delta that divides each dimension into layers, xi -> xi + delta
//      3) create arrays for each dimension, layerX,Y,Z ;
//              layerX has indices 0 to int(Xmax-Xmin), etc, for nX, nY, nZ elements 
//      4) each element of a layer is either nil, or a bit string for all the atoms, with 1 for each atom included in the layer
//      5) max storage is (nX + nY + nZ) * (nAtoms / 8 bytes)
//      6) take the midpoint of interval xmi ( = e.g. (x,i + x,i+1)/2.); position of any atom in interval is xmi +/- delta/2.,
//              distance dxij <= |xmi - xmj| + delta = |i - j|*delta + delta = (deltaNX + 1)*delta
//      7) so for a given 1D radius R, and any dimension, imin = ceil(i - R/delta + 1), imax = floor(i + R/delta - 1)
//      8) for simplicity, if we need simultaneous application of radius test in multiple dimensions, take union of objects in all directions
//              ^^ this overcounts, effectively using manhattan metric 
//              in 2D can do better, let deltaNX(max) = R/delta - 1
//                  then for each deltaNX in range -deltaNX(max) to +deltaNX(max), 
//                  deltaNY(max) = floor( sqrt((R/delta)**2 - (deltaNX + 1)**2) - 1 ), sqrt argument > 0 
//              Have not worked it out, but pretty sure in 3D,
//                  deltaNZ(max) = floor( sqrt((R/delta)**2 - (deltaNX + 1)**2 - (deltaNY + 1)**2 ) - 1 ), sqrt argument > 0
//               


struct objectBits {

    enum objectBitError: Error {
    case sizeError

    
    }

    var bits: [UInt8]
    
    init(_ n:Int, _ turnOn:[Int]? ) {
        let nbytes = Int(ceil( Double(n) / 8.0 ))
        bits = [UInt8](repeating:0, count:n )

        if turnOn != nil {
            for oidx in turnOn! {
                let byteidx = Int( oidx / 8 )
                let bitidx = oidx % 8 
                bits[byteidx] |= 1 << bitidx
            }
        }
    }

    mutating func addIndex(_ index:Int) {
        let byteidx = Int( index / 8 )
        let bitidx = index % 8
        bits[byteidx] |=  1 << bitidx
    }

    mutating func union(_ with:objectBits) throws {

        if with.bits.count != bits.count {
            throw objectBitError.sizeError
        }

        bits = zip(bits, with.bits).map { $0 | $1 }  
    }

    mutating func intersection(_ with:objectBits) throws  {

        if with.bits.count != bits.count {
            throw objectBitError.sizeError
        }

        bits = zip(bits, with.bits).map { $0 & $1 }  
    }

    func indices() -> [Int] {
        var ON = [Int]()
        for (nb,b) in bits.enumerated() {
            for p in 0..<8 {
                if b & 1 << p != 0 {
                    ON.append(8*nb + p)
                }
            }
        }

        return ON
    }

    mutating func setBits(_ content:[UInt8]) {
        bits = content
    }

}


enum AXES:Int {
    case X = 0
    case Y = 1
    case Z = 2
}


struct LAYERS {
    var axis:AXES

    var mincoord:Double
    var maxcoord:Double

    var delta:Double

    var objects:[Any]

    var layers:[Int]

    var layerBits:[objectBits?]


    // 
    init( _ objects:[Any], _ layers:[Int], _ delta:Double, _ axis:AXES, _ minCoord:Double ) {
        self.objects = objects
        self.layers = layer
        let maxLayer = layers.max()
        let nLayers = maxLayer + 1
        self.delta = delta      
        self.axis = axis 
        self.mincoord = minCoord 
        self.maxcoord = minCoord + maxLayer*delta

        let numobj = self.objects.count

        layerBits = Array(repeating:nil, count:nlayers)

        for i in 0..<numobj {
            let lidx = self.layers[i]
            if layerBits[lidx] == nil {
                layerBits[lidx] = objectBits( numobj, [i] )
            }
            else {
                layerBits[lidx]!.addIndex(i)
            }
        }

    }


}