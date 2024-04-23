
import Foundation
import MathTools

import Dispatch

let computeQueue = DispatchQueue( label:"compute", attributes: .concurrent )
let blocksQueue = DispatchQueue( label:"blocks" )


struct AtomCircle {

    var center:Vector
    var radius:Double
    var exposure = [exposedArc]()
    var axis:AXES
    var index
    var removed:Bool

    init(_ center:[Double], _ radius:Double, _ index:Int, _ axis:AXES ) {
        self.center = position 
        self.radius = radius
        self.axis = axis
        self.index = index
        removed = false
        arcs.append(exposedArc(nil, nil, self))
    }


    func intersectArcs(_ A:exposedArc, _ B:exposedArc ) -> exposedArc? {
        if A.ustart == nil {
            return B
        }

        if B.ustart == nil {
            return A
        }

        let AsInB = B.contains(A.ustart!)
        let AeInB = B.contains(A.uend!)
        if AsInB {
            if AeInB {
                return exposedArc(A.ustart!, A.uend!, self)
            }
            else {
                return exposedArc(A.ustart!, B.uend!, self)
            }
        }

        let BsInA = A.contains(B.ustart!)
        let BeInA = A.contains(B.uend!)

        if BsInA {
            if BeInA {
                return exposedArc(B.ustart!, B.uend!, self)
            }
            else {
                return exposedArc(B.ustart!, A.uend!, self)
            }
        }

        // no intersection

        return nil 


    }

    mutating func updateExposure(_ newarc:exposedArc) {
        // 
        var newExposure = [exposedArc]()

        for arc in exposure {
            let inter = intersectArcs(newarc, arc)
            if inter != nil {
                newExposure.append(inter!)
            }
        }

        self.exposure = newExposure

        if newExposure.count == 0 {
            self.removed = true
        }
    } 



    
}

struct Contour() {
    public init(_ circles:[AtomCircle]) {
        // assume all have the same axis
        let axis = circles[0].axis 

        
    }


}

// ref axis is just unit vector along X, Y or Z

let axisREF = [Vector(1.0,0.0,0.0), Vector(0.0,1.0,0.0), Vector(0.0,0.0,1.0)]

// 'right' and 'up' are references for planes perpendicular to X Y OR Z
// X -> Y,Z, Y -> Z,X , Z -> X,Y

let axisRIGHT = [1, 2, 0]
let axisUP    = [2, 0, 1]


struct exposedArc {

    // pstart, pend, etc = nil means that entire circle is exposed arc

    var acute:Bool 
    var dotSE:Double
    var pstart:Vector?
    var pend:Vector?

    var ustart:Vector?
    var uend:Vector?

    var atomcircle:AtomCircle


    init( _ start:Vector?, _ end:Vector?, _ atomcircle:AtomCircle) {
        if(start == nil) {
            pstart = nil
            ustart = nil
            pend = nil
            uend = nil
            acute = false
            dotSE = 0.0
            self.atomcircle = atomcircle
            
        }
        else {

            pstart = start
            pend = end 
            self.atomcircle = atomcircle

            ustart = pstart!.diff(self.atomcircle.center)
            ustart = ustart!.scale(1.0/ustart!.length())

            uend! = pend!.diff(self.atomcircle.center)
            uend! = uend!.scale(1.0/uend!.length())

            let ref = axisREF[self.atomcircle.axis.rawValue]

            acute = ref.dot(uend!.cross(with:ustart!)) > 0.0

            dotSE = ustart!.dot(uend!)
        }

    }

    func contains( _ u:Vector) -> Bool {
        if ustart == nil {
            return true 
        }

        let dotS = ustart!.dot(u)
        let dotE = uend!.dot(u)

        if min(dotS,dotE) > dotSE {
            if acute {
                return false 
            }
            else {
                return true
            }
        }
        else {
            if acute {
                return true 
            }
            else {
                return false
            }
        }
    }


    
}

// find intersection between circles, update exposed arcs and/or remove buried circle

func intersectAtomCircles(_ circleA:atomCircle, _ circleB:atomCircle) {
    let uAB = circleB.center.diff(circleA.center)
    let dAB = uAB.length()
    let rA = circleA.radius
    let rB = circleB.radius
    // no intersection 

    if dAb > (rA + rB) {
        return
    }

    // B in A 
    if dAB + circleB.radius < circleA.radius {
        circleB.remove = true 
        return
    }

    // A in B 
    if dAB + circleA.radius < circleB.radius {
        circleA.remove = true 
        return
    }

    let RIGHT = axisRIGHT[circleA.axis.rawValue]
    let UP = axisUP[circleA.axis.rawValue]
    uAB = uAB.scale(1.0/dAB)
    vAB = Vector([0.0,0.0,0.0])
    vAB[RIGHT] = -uAB[UP]
    vAB[UP] = uAB[RIGHT]

    let m = (dAB*dAB + rA*rA - rB*rB)/2.0

    let h = sqrt(rA*rA - m*m)

    let pdown = circleA.center.add(uAB.scale(m)).add(vAB.scale(-h))
    let pup = circleA.center.add(uAB.scale(m)).add(vAB.scale(h))

    // clockwise ; for A, pdown to pup, for B, pup to pdown

    circleA.updateExposure(exposedArc(pdown, pup, circleA))
    circleB.updateExposure(exposedArc(pup, pdown, circleB))


}

func circlesForAtoms( atompos:Matrix<Double>, radii:[Double], proberad:Double,
       axis:AXES, minCoord:Double, delta:Double, limits:[Int,Int], thread:Int ) 
        -> ([(AtomCircle,Int)],Int)   {
        
    var data = [(AtomCircle,Int)]()

    for aidx in limits[0]..<limits[1] {

        let catom = (0..<3) .map { atompos.storage[3*aidx + $0]}
        let caxis = catom[axis.rawValue]

        let augrad = radii[aidx] + proberad

        // bottom layer is below bottom of augmented sphere, so start one higher
        let layerlow = Int(((caxis - augrad) - minCoord)/delta) + 1
        // upper layer contains top of augmented sphere, include it in range
        let layerhi = Int(((caxis + augrad) - minCoord)/delta) + 1

        for layer in layerlo..<layerhi {
            layerc = minCoord + layer*delta 
            let h = abs(caxis - layerc) 
            // sanity test
            if h > augrad {
                continue
            }
            let crad = sqrt(augrad*augrad - h*h)
            data.append((AtomCircle(catom, crad, aidx, axis ),layer))
        }
    }

    return (data,thread)
        
}

func addBLOCK( _ BLOCKS: inout [[(AtomCircle,Int)]?], data:([(AtomCircle,Int)],Int) ) {
    BLOCK[data.1] = data.0
}

func atomCirclesForLayers( atompos:Matrix<Double>, radii:[Double], 
    proberad:Double, layerdelta:Double, axis:AXES, numthreads:Int=1 ) -> ([LAYER],[[AtomCircle]]) {

    let maxRad = radii.max()
    let touchRadius = maxRad + proberad
    let numatoms = atompos.shape[0]

    let axiscoords = (0..<numatoms) .map { atompos.storage[3*$0 + axis.rawValue] }

    let minaxiscoord = axiscoords.min() - proberad 

    var LIMITS = [[Int,Int]]()

    let atomsPerChunk = Int(atompos.shape[0]/numthreads)
    let nchunks = Int(ceil(Double(atompos.shape[0])/Double(atomsPerChunk)))

    for ichunk in 0..<nchunks {
        let lo = ichunk*atomsPerChunk
        let hi = lo + atomsPerChunk
        if ichunk == nchunks - 1 {
            hi = numatoms
        }

        LIMITS.append([lo,hi])
    }  

    // each block will hold list of atom circles coupled with layer index

    let BLOCKS:[[(AtomCircle,Int)]?] = Array(repeating:nil, count:nchunks)

    let group = DispatchGroup() 

    for ichunk in 0..<nchunks {

        group.enter() 

        computeQueue.async {
            let data = circlesForAtoms( atompos:atompos, radii:radii, proberad:proberad,
                axis:axis, minCoord:minaxiscoord, delta:delta, limits:LIMITS[ichunk], 
                thread:ichunk ) 

            blocksQueue.sync {
                addBlock( &BLOCKS, data )
            }

            group.leave()

        }

        

    }

    group.wait()

    var circles = [AtomCircle]()
    var circleLayers = [Int]() 

    for block in BLOCKS {
        circles += block! .map { b.0 }
        circlelayers += block! .map { b.1 }
    }

    let circleLAYERS = LAYERS(circles, circlelayers, delta, axis, minaxiscoord )

    return circleLAYERS
}

func findContours( _ circleLAYERS:LAYERS ) ->  {

    // get intersecting circles for each layer

    var radii = [Double]()
    var centercoord = [Double]()

    for o in circleLAYERS.objects {
        let circle = o as! AtomCircle 
        radii.append(circle.radius)
        centercoord += circle.center.coords
    }

    let centerCoord = Matrix<Double>([circleLAYERS.objects.count,3], content:centercoord)
    let Radii = Matrix<Double>([circleLAYERS.objects.count], content:radii)


    var dists:Matrix<Double>?

    do {
        try dists = cdist(centerCoord,centerCoord, numthreads:10)
        try dists.setdiagonal(1.0e12)
    }
    catch {
        print("unexpected exception in cdist/setdiagonal")
        return nil
    }








}

