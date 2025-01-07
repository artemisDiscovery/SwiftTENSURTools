
import Foundation
import MathTools
import SwiftMC33Lib

import Dispatch

let computeQueue = DispatchQueue( label:"compute", attributes: .concurrent )
let blocksQueue = DispatchQueue( label:"blocks" )



// try str.write(to: URL(fileURLWithPath:"probes.txt"), atomically: true, encoding: String.Encoding.utf8)

public func exportAtoms( _ path:String, _ atompos:[Vector], _ radii:[Double], _ rgb:[Vector]?=nil ) throws {
    var outstr = ""

    for (aidx,coord) in atompos.enumerated() {
        var r = 1.0
        var g = 1.0
        var b = 1.0
        if rgb != nil {
            r = rgb![aidx].coords[0]
            g = rgb![aidx].coords[1]
            b = rgb![aidx].coords[2]
        }
        outstr += "\(atompos[aidx].coords[0]) \(atompos[aidx].coords[1]) \(atompos[aidx].coords[2]) \(radii[aidx]) \(r) \(g) \(b)\n"
    }

    do {
        try outstr.write(to: URL(fileURLWithPath:path), atomically: true, encoding: String.Encoding.utf8)
    }
    catch {
        print("could not write spheres file \(path)")
    }
}





public class AtomCircle {

    var center:Vector
    var radius:Double
    var exposure = [exposedArc]()
    var axis:AXES
    var atom:Int
    var removed:Bool


    init(_ atom:Int, _ center:Vector, _ radius:Double,  _ axis:AXES ) {
        self.atom = atom 
        self.center = center 
        self.radius = radius
        self.axis = axis
    
        self.removed = false
  
        
    }

    func setRemoved(_ state:Bool) {
        self.removed = state
    }


    func updateExposure( _ pstart:Vector, _ pend:Vector, _ othercircle: AtomCircle) {


        let newarc = exposedArc(pstart, pend, othercircle, othercircle, self )
        

        if exposure.count == 0 && !removed {
            exposure.append(newarc)
            return 
        }
  
        var newexposure = [exposedArc]()

        for iarc in 0..<exposure.count {
            newexposure += intersectArcs( exposure[iarc], newarc )
        }

        exposure = newexposure

        if exposure.count == 0 {
            self.removed = true
        }


    }

    func intersectArcs(_ A: exposedArc, _ B: exposedArc ) -> [exposedArc] {
        
        let AstartInB = B.contains(A.ustart)
        let AendInB   = B.contains(A.uend)

        var returnArcs = [exposedArc]() 

        

        if AstartInB {
            if AendInB {
                // check for Astart, Aend in clockwise order in B
                if !B.clockwiseOrder(A.ustart, A.uend) {
                    let arc0 = exposedArc(B.pstart, A.pend, B.atomcircleStart, A.atomcircleEnd, self)
                    //arc0.setEndCircle(A.atomcircleEnd)
                    let arc1 = exposedArc(A.pstart, B.pend, A.atomcircleStart, B.atomcircleEnd, self)
                    //arc1.setEndCircle(B.atomcircleEnd)

                    returnArcs.append( arc0 )
                    returnArcs.append( arc1 )
                }
                else {
                    returnArcs.append( A )
                    
                }
                return returnArcs
            }
            else {
                // Astart in B and Aend NOT in B
                let arc0 = exposedArc(A.pstart, B.pend, A.atomcircleStart, B.atomcircleEnd, self)
                //arc0.setEndCircle(B.atomcircleEnd)
                returnArcs.append( arc0 )
                
                return returnArcs
            }
        }
        else if AendInB {
            // Astart NOT in B, Aend in B
            let arc0 = exposedArc(B.pstart, A.pend, B.atomcircleStart, A.atomcircleEnd, self)
            //arc0.setEndCircle(A.atomcircleEnd)
            returnArcs.append( arc0 )
            
            return returnArcs 
        }

        let BstartInA = A.contains(B.ustart)
        let BendInA   = A.contains(B.uend)

        if BstartInA {
            if BendInA {
                // check clockwise order Bstart, Bend
                if !A.clockwiseOrder( B.ustart, B.uend ) {
                    let arc0 = exposedArc(A.pstart, B.pend, A.atomcircleStart, B.atomcircleEnd, self)
                    //arc0.setEndCircle(B.atomcircleEnd)
                    let arc1 = exposedArc(B.pstart, A.pend, B.atomcircleStart, A.atomcircleEnd, self)
                    //arc1.setEndCircle(A.atomcircleEnd)
                    returnArcs.append( arc0 )
                    returnArcs.append( arc1 )

                }
                else {
                    returnArcs.append( B )
                }
                return returnArcs
            }
            else {
                // Bstart in A and Bend NOT in A
                let arc0 = exposedArc(B.pstart, A.pend, B.atomcircleStart, A.atomcircleEnd, self)
                //arc0.setEndCircle(A.atomcircleEnd)
                returnArcs.append( arc0 )

                return returnArcs
            }
        }
        else if BendInA {
            // Bstart NOT in A, Bend in A
            let arc0 = exposedArc(A.pstart, B.pend, A.atomcircleStart, B.atomcircleEnd, self)
            //arc0.setEndCircle(B.atomcircleEnd)
            returnArcs.append( arc0 )
            return returnArcs 
        }

        // no intersection

        return returnArcs 


    }

    
    public func str() -> String {
        var returnString = ""

        returnString += "center : \(center.coords[0]) , \(center.coords[1]), \(center.coords[2])\n"
        returnString += "radius : \(radius)\n"
        returnString += "axis : \(axis.rawValue)\n"
        returnString += "atom : \(atom)\n"
        returnString += "removed : \(removed)\n"
        returnString += "# exposed arcs : \(exposure.count)\n"

        for arc in exposure {
            returnString += arc.str()
        }

        return returnString
    }


    
}

public enum ContourError: Error {
    case missingCircleError
    case missingArcError
    case noInitialArc
    
    
}

let printArc = 
"""
def printArc( centerX , centerY, radius, color, ustartX, ustartY, uendX, uendY, ax, num=10 ) :
    # clockwise arc, but angles handle CCW, make uend the 'start' of the drawn arc
    center = np.array([centerX,centerY])
    ref = np.array([uendX,uendY])
    perp = np.array([0.0,0.0])
    perp[0] = -ref[1]
    perp[1] = ref[0]
    ustart = np.array([ustartX,ustartY])
    # frame coords of ustart
    ustart_ref = np.dot(ustart,ref) 
    ustart_perp = np.dot(ustart,perp)
    #
    angle = math.atan2( ustart_perp, ustart_ref )
    if angle < 0. : angle += 2.*math.pi
    delta = angle / num 
    #
    X = []
    Y = []
    #
    for i in range(num) :
        astart = i * delta
        aend = astart + delta 
        pstart = center + radius*math.cos(astart)*ref + radius*math.sin(astart)*perp
        pend = center + radius*math.cos(aend)*ref + radius*math.sin(aend)*perp
        if i == 0 :
            ax.plot([pstart[0],pend[0]],[pstart[1],pend[1]],color='magenta')
            continue
        X.append(pstart[0])
        X.append(pend[0])
        Y.append(pstart[1])
        Y.append(pend[1])
    ax.plot(X,Y,color=color)

"""

public func printPython( _ circles:[AtomCircle], _ arcs:[exposedArc]) -> String {
    // print out text to copy into python session 
    //
    // Show removed circles in yellow, retained circles in red, retained singleton in black,
    // contour arcs in blue, non-contour arcs in green

    var text = "import matplotlib.pyplot as plt\nplt.ion()\nimport numpy as np\nimport math\n"
    text += printArc
    text += "\n"

    // get overall dimensions 

    let axis = circles[0].axis
    let RIGHT = axisRIGHT[axis.rawValue]
    let UP = axisUP[axis.rawValue]

    // use  RIGHT and UP as 'x' and 'y' 

    var xmin = circles .map { $0.center.coords[RIGHT] - $0.radius } .min()!
    xmin -= 1.0
    var xmax = circles .map { $0.center.coords[RIGHT] + $0.radius } .max()!
    xmax += 1.0
    var ymin = circles .map { $0.center.coords[UP] - $0.radius } .min()!
    ymin -= 1.0
    var ymax = circles .map { $0.center.coords[UP] + $0.radius } .max()!
    ymax += 1.0

    let deltaX = xmax - xmin 
    let deltaY = ymax - ymin 
    let xmid = (xmax + xmin)/2.0
    let ymid = (ymax + ymin)/2.0

    var XMIN = 0.0
    var XMAX = 0.0
    var YMIN = 0.0
    var YMAX = 0.0 

    if deltaX > deltaY {
        XMIN = xmin
        XMAX = xmax 
        YMIN = ymid - (deltaX/2.0)
        YMAX = ymid + (deltaX/2.0)
    }
    else {
        YMIN = ymin
        YMAX = ymax 
        XMIN = xmid - (deltaY/2.0)
        XMAX = xmid + (deltaY/2.0)
    }

    text += "fig, ax = plt.subplots()\n"
    text += "ax.set_xlim(\(XMIN), \(XMAX))\n"
    text += "ax.set_ylim(\(YMIN), \(YMAX))\n"
    text += "ax.set_aspect('equal')\n"

    var color = ""

    for circle in circles {
        if circle.removed {
            color = "yellow"
        }
        else if circle.exposure.count == 0 {
            color = "black"
        }
        else {
            color = "red"
        }

        text += "p = plt.Circle((\(circle.center.coords[RIGHT]),\(circle.center.coords[UP])), \(circle.radius), color='\(color)', fill=False)\n"
        text += "ax.add_patch(p)\n"

        if !circle.removed {
            for arc in circle.exposure {
            text += "printArc( \(circle.center.coords[RIGHT]) , \(circle.center.coords[UP]), \(circle.radius), 'green', \(arc.ustart.coords[RIGHT]), \(arc.ustart.coords[UP]), \(arc.uend.coords[RIGHT]), \(arc.uend.coords[UP]), ax, num=10 )\n"
            }
        }
        
    }


    for arc in arcs {
        let circle = arc.parentcircle
        text += "printArc( \(circle.center.coords[RIGHT]) , \(circle.center.coords[UP]), \(circle.radius), 'blue', \(arc.ustart.coords[RIGHT]), \(arc.ustart.coords[UP]), \(arc.uend.coords[RIGHT]), \(arc.uend.coords[UP]), ax, num=10 )\n"
    }

    return text

}


public struct Contour {

    var axis:AXES
    var arcsInOrder = [exposedArc]() 
    var singletons = [AtomCircle]()
    var circles = [AtomCircle]()
    var clockwise:Bool

    public init(_ circles: [AtomCircle]) throws  {

        self.circles = circles
        // assume all have the same axis
        axis = circles[0].axis 

        var initialArc:exposedArc? = nil 
        

        for circle in circles {
            if circle.removed  {
                continue
            }

            if circle.exposure.count == 0 {
                singletons.append(circle)
                circle.setRemoved(true)
            }
        }

        for cidx in 0..<circles.count {
            let circle = circles[cidx]
                    
            if circle.removed  {
                continue 
            }


            for arc in circle.exposure {
                if arc.removed {
                    continue
                }

                initialArc = arc 
                
                break
            }

            if initialArc != nil {
                break
            }

        }

        if initialArc == nil {

            
            //print("contour: no initial arc")
            // 
            if singletons.count == 0 {
                throw ContourError.noInitialArc
            }
            
        }


        var currentArc:exposedArc?
        

        //print("initial circle: ")
        //print(initialCircle!.str())
        //print("initial arc: ")
        //print(initialArc!.str())
        if initialArc != nil {
            arcsInOrder.append(initialArc!)
            currentArc = initialArc
        }
        
        
        var closed = false
        var accumAngle = 0.0
        var prevdisp:Vector? = nil

        if initialArc == nil  {
            closed = true
        }

        while !closed {
        //for _ in 0..<10000 {
            //print("contour: current arc :")
            //print(currentArc.str())
            let nextCircle = currentArc!.atomcircleEnd
            
            if nextCircle.removed {
                print("contour: expected circle removed")
                //print("current circle :")
                //print(currentArc!.parentcircle.str())
                //print("current arc :")
                //print(currentArc!.str())
                //print("missing circle :")
                //print(nextCircle.str())
                throw ContourError.missingCircleError
            }

            var bestArc:exposedArc? = nil
            var bestDist:Double = Double.infinity

            //print("contour: next circle :")
            //print(nextCircle.str())

            for arc in nextCircle.exposure {
                //print("contour: test arc :")
                //print(arc.str())
                if arc.removed {
                    continue
                }
                let d = arc.pstart.sub(currentArc!.pend).length()

                if d < bestDist {
                    bestArc = arc 
                    bestDist = d
                }

            }

            if bestArc == nil {
                print("contour: no best next arc found")
                throw ContourError.missingArcError
            }


            arcsInOrder.append(bestArc!)

           
            if bestArc!.pend == initialArc!.pstart {
                closed = true
            }
        

            let disp = bestArc!.pstart.sub(currentArc!.pstart).unit()

            if prevdisp != nil {
                var ang = acos(disp!.dot(prevdisp!))
                // define positive angle as clockwise
                // cross product expensive, better way?
                if disp!.cross(prevdisp!).dot(axisREF[axis.rawValue]) < 0.0 {
                    ang = -ang
                }

                accumAngle += ang
            }

            currentArc = bestArc!
            prevdisp = disp 
        }

        // assume closed here
        clockwise = true 

        if arcsInOrder.count > 2 {
            clockwise = accumAngle > 0.0
        }
        

        for aidx in 0..<self.arcsInOrder.count {
            arcsInOrder[aidx].setRemoved(true)
        }

        
    }

    public func str(verbose:Bool=false) -> String {

        var report = ""

        report += "# singletons = \(singletons.count)\n"

        if verbose {
            for circle in singletons {
                report += circle.str()
            }
        }
        else {
            for circle in singletons {
                report += "\tatom \(circle.atom)\n"
            }
        }

        report += "\n# arcs = \(arcsInOrder.count)\n"

        if verbose {
            for arc in arcsInOrder {
                report += arc.str()
            }
        }
        else {
            for arc in arcsInOrder {
                report += "\tstart atom : \(arc.atomcircleStart.atom), end atom : \(arc.atomcircleEnd.atom)\n"
            }
        }

        report += "\nclockwise : \(clockwise)"

        return report 

    }


}

// ref axis is just unit vector along X, Y or Z

let axisREF = [Vector([1.0,0.0,0.0]), Vector([0.0,1.0,0.0]), Vector([0.0,0.0,1.0])]

// 'right' and 'up' are references for planes perpendicular to X Y OR Z
// X -> Y,Z, Y -> Z,X , Z -> X,Y

let axisRIGHT = [1, 2, 0]
let axisUP    = [2, 0, 1]

let translateAXIS = [ AXES.X, AXES.Y, AXES.Z ]


public class exposedArc {

    // pstart, pend points to define

    var acute:Bool 
    var dotSE:Double
    var clockwiseMeasure:Double
    var pstart:Vector
    var pend:Vector

    var ustart:Vector
    var uend:Vector


    var parentcircle:AtomCircle
    

    // intersections labeled by intersecting circle  
    var atomcircleStart:AtomCircle
    var atomcircleEnd:AtomCircle

    var removed:Bool
  

    
    init( _ start:Vector, _ end:Vector, _ startcircle: AtomCircle, _ endcircle: AtomCircle, _ parentcircle:AtomCircle ) {
        
            pstart = start
            pend = end 
            atomcircleStart = startcircle
            atomcircleEnd = endcircle

            self.parentcircle = parentcircle

            //self.axis = parentcircle.axis
            //self.center = parentcircle.center
            


            ustart = start.sub(parentcircle.center)
            ustart = ustart.scale(1.0/ustart.length())

            uend = end.sub(parentcircle.center)
            uend = uend.scale(1.0/uend.length())

            let ref = axisREF[parentcircle.axis.rawValue]

            acute = ref.dot(uend.cross(ustart)) > 0.0

            dotSE = ustart.dot(uend)

            // I can't call the function ??

            // clockwiseMeasure = clockwise(uend)

            let crossdot = uend.cross(ustart).dot(ref)

            if crossdot >= 0 {
                clockwiseMeasure =  1.0 - uend.dot(ustart)
            }
            else {
                clockwiseMeasure =  3.0 + uend.dot(ustart)
            }

                removed = false 

 

    }

    func setEndCircle( _ end: AtomCircle) {
        self.atomcircleEnd = end
    }

    func setRemoved(_ state:Bool) {
        self.removed = state 
    }


    // return inside/outside for two vectors

    func contains( _ u:Vector) -> Bool {
        

        let dotS = ustart.dot(u)
        let dotE = uend.dot(u)

        // cosine monotonic - if max is greater than dotSE, then vector is inside smaller angle
        // BUT, this fast test ONLY WORKS for 'smaller angle' < 120 degrees; otherwise vector at (360 - theta)/2 can violate test

        // 
        if dotSE > -0.5 {

            if min(dotS,dotE) > dotSE {
                // inside smaller angle - note that 'acute' is not properly used here, I meant '< 180 deg'
                if acute {
                    return true
                }
                else {
                    return false
                }
            }
            else {
                // outside smaller angle
                if acute {
                    return false 
                }
                else {
                    return true
                }
            }
        }
        
        // need to do slow test

        return clockwise(u) < clockwiseMeasure

        
    }

    func clockwise( _ u:Vector  ) -> Double {
        // return positive measure of clockwise angle w.r.t start of arc

        let ref = axisREF[parentcircle.axis.rawValue]

        let crossdot = u.cross(ustart).dot(ref)

        if crossdot >= 0 {
            return 1.0 - u.dot(ustart)
        }
        else {
            return 3.0 + u.dot(ustart)
        }

    }

    func clockwiseOrder( _ uA:Vector, _ uB:Vector ) -> Bool {

        let measureA = clockwise( uA )
        let measureB = clockwise( uB )

        return measureA < measureB 
    }

    public func str() -> String {
        var returnString = ""

        returnString += "\tarc pstart = \(pstart.coords[0]), \(pstart.coords[1]), \(pstart.coords[2])\n"
        returnString += "\tarc pend   = \(pend.coords[0]), \(pend.coords[1]), \(pend.coords[2])\n"
        returnString += "\tarc ustart = \(ustart.coords[0]), \(ustart.coords[1]), \(ustart.coords[2])\n"
        returnString += "\tarc uend   = \(uend.coords[0]), \(uend.coords[1]), \(uend.coords[2])\n"
        returnString += "\tparent circle atom = \(parentcircle.atom)\n"
        returnString += "\tstart intersection atom = \(atomcircleStart.atom)\n"
        returnString += "\tend intersection atom = \(atomcircleEnd.atom)\n"
        returnString += "\tdotSE : \(dotSE)\n"
        returnString += "\tacute : \(acute)\n"
        returnString += "\tremoved : \(removed)\n"
        


        return returnString
    }
    
}

// find intersection between circles, update exposed arcs and/or remove buried circle

public func intersectAtomCircles(_ circleA: AtomCircle, _ circleB: AtomCircle) {

    if circleA.removed && circleB.removed {
        return
    }


    var uAB = circleB.center.sub(circleA.center)
    let dAB = uAB.length()
    let rA = circleA.radius
    let rB = circleB.radius
    // no intersection 

    if dAB > (rA + rB) {
        //print("no intersection")
        return
    }

    // B in A 
    if dAB + circleB.radius < circleA.radius {
        circleB.setRemoved(true) 
        //print("remove circle \(circleB.atom), inside \(circleA.atom)")
        return
    }

    // A in B 
    if dAB + circleA.radius < circleB.radius {
        circleA.setRemoved(true)
        //print("remove circle \(circleA.atom), inside \(circleB.atom)")
        return
    }

    let RIGHT = axisRIGHT[circleA.axis.rawValue]
    let UP = axisUP[circleA.axis.rawValue]
    uAB = uAB.scale(1.0/dAB)
    var vAB = Vector([0.0,0.0,0.0])
    vAB.coords[RIGHT] = -uAB.coords[UP]
    vAB.coords[UP] = uAB.coords[RIGHT]

    let m = (dAB*dAB + rA*rA - rB*rB)/(2.0 * dAB)

    // check for valid intersection

    let disc = rA*rA - m*m 

    if disc <= 0.0 {
        return 
    }

    let h = sqrt(disc)

    let pdown = circleA.center.add(uAB.scale(m)).add(vAB.scale(-h))
    let pup = circleA.center.add(uAB.scale(m)).add(vAB.scale(h))

    // clockwise ; for A, pdown to pup, for B, pup to pdown

    circleA.updateExposure(pdown, pup, circleB)
    circleB.updateExposure(pup, pdown, circleA)


}



public func circlesForAtoms( atompos:Matrix<Double>, radii:[Double], proberad:Double,
       axis:AXES, minCoord:Double, delta:Double, limits:[Int], thread:Int ) 
        -> ([(AtomCircle,Int)],Int)   {
        
    var data = [(AtomCircle,Int)]()
    let storage = atompos.getStorage()
    let RIGHT = axisRIGHT[axis.rawValue]
    let UP = axisUP[axis.rawValue]

    var units:[Vector]? = nil
    var sizes:[Double]? = nil


    for aidx in limits[0]..<limits[1] {
        
        let catom = Vector((0..<3) .map { storage[3*aidx + $0]})
        let caxis = catom.coords[axis.rawValue]

        let augrad = radii[aidx] + proberad

        // bottom layer is below bottom of augmented sphere, so start one higher
        let layerlo = Int(( (caxis - augrad) - minCoord) / delta) + 1
        // upper layer contains top of augmented sphere, include it in range
        let layerhi = Int(((caxis + augrad) - minCoord)/delta) + 1


        for layer in layerlo..<layerhi {
            let layerc = minCoord + Double(layer)*delta 
            let h = abs(caxis - layerc) 
            // sanity test
            if h > augrad {
                continue
            }
            let crad = sqrt(augrad*augrad - h*h)

            var center = Vector([0.0,0.0,0.0])
            center.coords[axis.rawValue] = layerc
            center.coords[RIGHT] = catom.coords[RIGHT]
            center.coords[UP] = catom.coords[UP]


            // all atoms are in unit cell

            let theCircle = AtomCircle(aidx, center, crad, axis )
            data.append((theCircle,layer))


        }
    }


    return (data,thread)
        
}

func addBLOCK( _ BLOCKS: inout [[(AtomCircle,Int)]?], _ data:([(AtomCircle,Int)],Int) ) {
    BLOCKS[data.1] = data.0
}

public func atomCirclesForLayers( atompos:Matrix<Double>, radii:[Double], 
    proberad:Double, minaxiscoord:Double, layerdelta:Double, axis:AXES, numthreads:Int=1) -> LAYERS {

    
    let shape = atompos.getShape()
    let numatoms = shape[0]
    

    var LIMITS = [[Int]]()

    let atomsPerChunk = Int(shape[0]/numthreads)
    let nchunks = Int(ceil(Double(shape[0])/Double(atomsPerChunk)))

    for ichunk in 0..<nchunks {
        let lo = ichunk*atomsPerChunk
        var hi = lo + atomsPerChunk
        if ichunk == nchunks - 1 {
            hi = numatoms
        }

        LIMITS.append([lo,hi])
    }  

    // each block will hold list of atom circles coupled with layer index

    var BLOCKS:[[(AtomCircle,Int)]?] = Array(repeating:nil, count:nchunks)

    let group = DispatchGroup() 

    for ichunk in 0..<nchunks {

        group.enter() 

        computeQueue.async {
            let data = circlesForAtoms( atompos:atompos, radii:radii, proberad:proberad,
                axis:axis, minCoord:minaxiscoord, delta:layerdelta, limits:LIMITS[ichunk], 
                thread:ichunk ) 

            blocksQueue.sync {
                addBLOCK( &BLOCKS, data )
            }

            group.leave()

        }

        

    }

    group.wait()

    var circles = [AtomCircle]()
    var circleLayers = [Int]() 

    for block in BLOCKS {
        circles += block! .map { $0.0 }
        circleLayers += block! .map { $0.1 }
    }

    let circleLAYERS = LAYERS(circles, circleLayers, layerdelta, axis, minaxiscoord )

    return circleLAYERS
}

public func intersectCirclesInLayerRange( _ circleLayers:LAYERS, _ limits:[Int]) -> ([[AtomCircle]],[Int])? {

    let layerBits = circleLayers.layerBits

    let allcircles = circleLayers.objects as! [AtomCircle]

    var returnCircles = [[AtomCircle]]()
    var returnLayerIndices = [Int]()

    for lidx in limits[0]..<limits[1] {
        if layerBits[lidx] == nil {
            continue
        }

        let layercircles = layerBits[lidx]!.indices() .map { allcircles[$0] }

        if layercircles.count == 1 {
            returnCircles.append(layercircles)
            returnLayerIndices.append(lidx)
            continue
        }

        var centers = [Double]()
        var radii = [Double]()

        for circle in layercircles {
            centers += circle.center.coords
            radii.append(circle.radius)
        }

        let circleCenters = Matrix<Double>([layercircles.count,3], content:centers )
        let circleRadii = Matrix<Double>([layercircles.count], content:radii )

        var pairs:[[Int]]?

        do {
            let dists = try cdist(circleCenters,circleCenters)
            try dists.setdiagonal(1.0e12)
            let sumRadii = try circleRadii.addTranspose(circleRadii)

            let mask = try Mask.compare(dists, sumRadii) { $0 < $1 }
            pairs = mask.nonzero()
        }
        catch {
            print("exception in intersectCircles for layer \(lidx)")
            return nil 
        }


        for p in pairs! {
            if p[0] < p[1] {
                intersectAtomCircles( layercircles[p[0]], layercircles[p[1]] )
            }
        }

        returnCircles.append(layercircles)
        returnLayerIndices.append(lidx)


    }

    return (returnCircles,returnLayerIndices)

}


public func intersectingCirclesForLayers( _ circleLayers:LAYERS, probeRadius:Double, numthreads:Int=1, skipCCWContours:Bool=false ) -> ([[Contour]],[Probe]) {

    // set limits using fraction of total number of circles

    var LIMITS = [[Int]]()

    let atomcircles = circleLayers.objects as! [AtomCircle]
    let layerBits = circleLayers.layerBits

    let circlesPerChunk = Int(atomcircles.count/numthreads)

    // lo and hi specify layer ranges

    var currentLo = 0
    var currentCircleCount = 0
     
    while currentLo < layerBits.count  {

        
        if layerBits[currentLo] != nil {
            currentCircleCount = layerBits[currentLo]!.indices().count
        }


        if currentLo == layerBits.count - 1 {
            break 
        }

        var lastHi = -1 

        for currentHi in (currentLo+1)..<layerBits.count {
            if layerBits[currentHi] != nil {
                currentCircleCount += layerBits[currentHi]!.indices().count
            }
            if currentCircleCount >= circlesPerChunk {
                LIMITS.append([currentLo,currentHi+1])
                currentCircleCount = 0 
                lastHi = currentHi
                break
            }
        }

        if lastHi < 0 {
            break
        }

        currentLo = lastHi + 1

    }

    // add last if needed 

    if currentCircleCount > 0 {
        LIMITS.append([currentLo,layerBits.count])
    }

    var cumcount = 0 

    for lim in LIMITS {
        //print("limit \(lim)")
        var ccount = 0 
        for lidx in lim[0]..<lim[1] {
            if layerBits[lidx] == nil {
                continue 
            }
            ccount += layerBits[lidx]!.indices().count
        }
        //print("\tcount for range = \(ccount)")
        cumcount += ccount
    }

    //print("\ncumulative count = \(cumcount), compare to \(atomcircles.count) circles expected")

    var BLOCKS = [([[AtomCircle]],[Int])]()

    let group = DispatchGroup() 

    for ichunk in 0..<LIMITS.count {

        group.enter() 

        computeQueue.async  {
            let data = intersectCirclesInLayerRange(circleLayers,LIMITS[ichunk]) 

            blocksQueue.sync {
                BLOCKS.append(data!)
            }

            group.leave()

        }

        

    }

    group.wait()

    BLOCKS = BLOCKS.sorted { $0.1[0] < $1.1[0] }

    var sortedLayers = [Int]()
    var sortedCircles = [[AtomCircle]]()

    for block in BLOCKS {
        sortedLayers += block.1
        sortedCircles += block.0
    }

    //for lidx in 0..<sortedLayers.count {
    //    print("sorted layer \(sortedLayers[lidx]), \(sortedCircles[lidx].count) circles")
    //    
    //}

    var contoursForLevels = [[Contour]]()

    var probes = [Probe]()

    for lidx in 0..<sortedLayers.count {
        var contours = [Contour]()
        //print("layer \(lidx) contours ...")
        while true {
            var cont:Contour?
            do {
                cont = try Contour(sortedCircles[lidx])
            }
            catch {
                if (error as! ContourError) == ContourError.noInitialArc {
                    break
                }
                print("exception in Contour : \(error)")
                break
            }
            
            if cont == nil {
                print("nil contour, continue")
                continue
            }
            if cont!.singletons.count == 0 && cont!.arcsInOrder.count == 0 {
                break
            }
            

            if skipCCWContours && !cont!.clockwise {
                continue
            }

            contours.append(cont!)
            probes += probesForContour(cont!, probeRadius:probeRadius)

        }

        contoursForLevels.append(contours)
        //print("layer \(lidx) \(contours.count) contours ...")
        var nclock = 0
        var ncounterclock = 0
        for contour in contours {
            if contour.clockwise {
                nclock += 1

            }
            else {
                ncounterclock += 1
            }
        }
        //print("\t\(nclock) clockwise, \(ncounterclock) counterclockwise")
        //for cont in contours {
        //    print(cont.str())
        //}
    }

    //let plidx = 23
    //print("contour \(plidx) ")

    //print(printPython(contoursForLevels[plidx][1].circles,contoursForLevels[plidx][1].arcsInOrder ))

    print("\nfound \(probes.count) probes")

    return (contoursForLevels,probes)

}

public class Probe {
    public var center:Vector
    public var proberadius:Double
    public var atoms:[Int]
    public var singleton:Bool
    public var clockwisecontour:Bool

    init(center:Vector, radius:Double, atoms:[Int], singleton:Bool, clockwise:Bool) {

        self.center = center 
        self.proberadius = radius
        self.atoms = atoms
        self.singleton = singleton
        self.clockwisecontour = clockwise
    }
}

public func probesForContour( _ contour:Contour, probeRadius:Double, minOverlap:Double=0.5 ) -> [Probe] {

    let RIGHT = axisRIGHT[contour.axis.rawValue]
    let UP = axisUP[contour.axis.rawValue]

    let minsep = probeRadius - (minOverlap * probeRadius)


    var probes = [Probe]()

    for arc in contour.arcsInOrder {


        let atomC = arc.parentcircle.atom
        let atomS = arc.atomcircleStart.atom 
        
        let ref = arc.uend 
        var perp = Vector([0.0,0.0,0.0])

        perp.coords[RIGHT] = -ref.coords[UP]
        perp.coords[UP] = ref.coords[RIGHT]

        let x = arc.ustart.dot(ref)
        let y = arc.ustart.dot(perp)
        
        var angle = atan2(y, x)

        if angle < 0.0 {
            angle += 2.0 * acos(-1.0)
        }

        let radius = arc.parentcircle.radius

        var num = max(Int(ceil((angle*radius)/minsep)),2)

        if (minsep/2.0) < radius {
            let mindelta = 2.0 * asin((minsep/2.0)/radius)

            num = max( Int(ceil(angle/mindelta)), 2)
        }
        

        let delta = angle/Double(num)

        

        // go to one position less than end, avoid duplicate probes 

        for j in 0..<num {
            var atoms = [atomC]
            let ang = angle - Double(j)*delta

            if j == 0 {
                atoms.append(atomS)
            }
            
            let pos = arc.parentcircle.center + (radius*cos(ang))*ref + (radius*sin(ang))*perp

            probes.append(Probe(center:pos, radius:probeRadius, atoms:atoms, singleton:false, clockwise:contour.clockwise))


        }
    }

    for circle in contour.singletons {

        let atoms = [circle.atom]
        let radius = circle.radius
        let num = max( Int((2.0*acos(-1.0)*radius)/minsep), 4 )
        let delta = (2.0*acos(-1.0))/Double(num)

        var ref = Vector([0.0,0.0,0.0])

        ref.coords[RIGHT] = 1.0 

        var perp = Vector([0.0,0.0,0.0])

        perp.coords[UP] = 1.0 

        // start at 2*PI - delta, go all way to zero

        // why did I do this??
        //let startang = 2.0*acos(-1.0) - delta

        for j in 0..<(num+1) {
            //let ang = startang - Double(j)*delta
            let ang = Double(j)*delta
            let pos = circle.center + (radius*cos(ang))*ref + (radius*sin(ang))*perp
            probes.append(Probe(center:pos, radius:probeRadius, atoms:atoms, singleton:true, clockwise:true))
        }

        

    }


    return probes 


}

// returns probes for all axes used, and levels for X, Y and Z

public func generateSurfaceProbes( coordinates:[Vector], radii:[Double], probeRadius:Double, levelspacing:Double, minoverlap:Double, numthreads:Int,
        skipCCWContours:Bool,  unitcell:UnitCell?=nil, atomindices:[Int]?=nil, debugAXES:[AXES]?=nil ) 
        -> ([Probe],[[Double]]) {


    let radiiVec = Vector(radii)

    var atomcoords = [Double]() 

    for vec in coordinates {
        atomcoords += vec.coords 
    }

    let atompos = Matrix<Double>([coordinates.count,3], content:atomcoords )

    var probes = [Probe]()
    var layers = [[Double]]()

    var useAXES = [ AXES.X, AXES.Y, AXES.Z ]

    if debugAXES != nil {
        useAXES = debugAXES!
    }

    
    if unitcell != nil {
        useAXES = []
        for ax in [ AXES.X, AXES.Y, AXES.Z ] {
            if ax != unitcell!.membraneaxis {
                useAXES.append(ax)
            }
        }
    }
    

    for (naxis,axis) in useAXES.enumerated() {

        print("Axis : \(axis.rawValue)")

        let axiscoords = Vector( coordinates .map { $0.coords[axis.rawValue] } )

        let lowermin = axiscoords - radiiVec
        
        
        let lowlim = lowermin.coords.min()! - probeRadius 

        var minaxiscoord = ceil(abs(lowlim)/levelspacing) * levelspacing 

        if lowlim < 0 {
            minaxiscoord = -minaxiscoord
        }

        let atomcircleLAYERS = atomCirclesForLayers( atompos:atompos, radii:radii, 
                proberad:probeRadius, minaxiscoord:minaxiscoord, layerdelta:levelspacing, axis:axis, numthreads:numthreads )
        
        let probedata = intersectingCirclesForLayers(atomcircleLAYERS, probeRadius:probeRadius, numthreads:numthreads, skipCCWContours:skipCCWContours )

        var axisprobes = probedata.1

        // if atom indices defined, translate probe atom indices to original indices

        if atomindices != nil {
            for probe in axisprobes { 
                probe.atoms = probe.atoms .map { atomindices![$0] }
            }
        }

        probes += axisprobes

        layers.append([atomcircleLAYERS.mincoord, atomcircleLAYERS.maxcoord])



    }

    return (probes,layers)
}

public enum TriangulationError: Error {
    case gridProbeDistanceError
    case densityError
    
}

public func indexFromIndices( _ indices:[Int],  shape:[Int], strides:[Int] )  -> Int {
    
    if indices.count != shape.count {
        return -1 
    }

    var index = 0 

    // .enumerated seems very slow

    for sidx in 0..<indices.count {
        let idx = indices[sidx]
        if idx < 0 || idx >= shape[sidx] {
            return -1
        }

        index += strides[sidx] * idx
    }

    return index 
}


// This returns for argument probe position a list of density values and 
// offsets in the original grid storage

public func densityForProbes( probes:[Probe], radius:Double, delta:Double, epsilon:Double, 
    griddeltas:[Double], limits:[[Double]], gridShape:[Int], gridStrides:[Int] ) -> [Double] {

    
        var linearDensity = Array( repeating:0.0, count:gridShape[0]*gridShape[1]*gridShape[2] )
        
        let gradii = griddeltas .map { Int(radius/$0) + 1 }

        let A = -1.0 / (delta + epsilon)
        let B = (radius + epsilon)/(delta + epsilon)

        for probe in probes {

            let center = probe.center

            let gcenter = (0..<3) .map { Int((center.coords[$0] - limits[$0][0])/griddeltas[$0]) }

            let gmin = (0..<3) .map { max(gcenter[$0] - gradii[$0], 0) }
            let gmax = (0..<3) .map { min(gcenter[$0] + gradii[$0], gridShape[$0]) }

            // get ranges, check for probe off grid

            for k in 0..<3 {
                if gmax[k] < gmin[k] {
                    continue
                }
            }

            var probeLinearCoords = [Vector]()
            var probeLinearIndices = [Int]()
            
            
            //print("gmin, gmax = \(gmin) , \(gmax)")

            for iz in gmin[2]...gmax[2] {
                let z = limits[2][0] + Double(iz)*griddeltas[2]
                for iy in gmin[1]...gmax[1] {
                    let y = limits[1][0] + Double(iy)*griddeltas[1]
                    for ix in gmin[0]...gmax[0] {
                        let offset = ix*gridStrides[0] + iy*gridStrides[1] + iz
                        
                        let x = limits[0][0] + Double(ix)*griddeltas[0]
                        probeLinearCoords.append(Vector([x,y,z]))
                        probeLinearIndices.append(offset)
                    }
                }

            }

            let dists = probeLinearCoords .map { $0.dist(probe.center) }

            

            for (d,idx) in zip(dists,probeLinearIndices) {
                
                if d < radius - delta {
                    linearDensity[idx] += 1.0
                    if linearDensity[idx] >= 1.0 {
                        linearDensity[idx] = 1.0
                    }
                }
                else if d <= radius + epsilon {
                    linearDensity[idx] += A*d + B
                    if linearDensity[idx] >= 1.0 {
                        linearDensity[idx] = 1.0
                    }
                }
            }



        }

        // 

        return linearDensity


    }


public func runMarchingCubes( density:Matrix<Double>, limits:[[Double]], griddeltas:[Double], gridvertices:[Int], isoLevel:Double ) 
        ->  ([Vector], [Vector], [[Int]]) {
    let nx:UInt32 = UInt32(gridvertices[0])
    let ny:UInt32 = UInt32(gridvertices[1])
    let nz:UInt32 = UInt32(gridvertices[2])

    //print("in runMarchingCubes, nx,ny,nz = \(nx), \(ny), \(nz)")

    let ptr_r0 = UnsafeMutablePointer<Double>.allocate(capacity: 3)
    ptr_r0.initialize(repeating: 0.0, count: 3)
        defer {
            ptr_r0.deinitialize(count: 3)
            ptr_r0.deallocate()
        }

    let r0 = UnsafeMutableBufferPointer(start:ptr_r0, count:3 )

    r0[0] = limits[0][0]
    r0[1] = limits[1][0]
    r0[2] = limits[2][0]

    //print("grid lower limits : \(r0[0]), \(r0[1]), \(r0[2])")

    let ptr_d = UnsafeMutablePointer<Double>.allocate(capacity: 3)
        ptr_d.initialize(repeating: 0.0, count: 3)
        defer {
            ptr_d.deinitialize(count: 3)
            ptr_d.deallocate()
        }

    let d = UnsafeMutableBufferPointer(start:ptr_d, count:3 )

    d[0] = griddeltas[0]
    d[1] = griddeltas[1]
    d[2] = griddeltas[2]
    
    //print("grid deltas : \(d[0]), \(d[1]), \(d[2])")
    
    let ptr = UnsafeMutablePointer<GRD_data_type>.allocate(capacity: density.storage.count)

    ptr.initialize(repeating: 0.0, count: density.storage.count)
    defer {
        ptr.deinitialize(count: density.storage.count)
        ptr.deallocate()
    }

    let data = UnsafeMutableBufferPointer(start:ptr, count:density.storage.count )

    _ = (0..<density.storage.count) .map { data[$0] = density.storage[$0]}

    let G:UnsafeMutablePointer<_GRD> = grid_from_data_pointer(nx, ny, nz, ptr_r0, ptr_d, ptr)

    let M:UnsafeMutablePointer<MC33> = create_MC33(G)

    let S:UnsafeMutablePointer<surface> = calculate_isosurface(M, Float(isoLevel))

    let nV = S.pointee.nV
    let nT = S.pointee.nT 

    //print("have \(nV) vertices, \(nT) faces")

    let V = S.pointee.V
    let N = S.pointee.N 
    let T = S.pointee.T


    var VERTICES = [Vector]()
    var NORMALS = [Vector]()

    // it appears that vertices are in grid coordinates !!

    for iv in 0..<Int(nV) {
        var x = Double(V![iv].0)
        var y = Double(V![iv].1)
        var z = Double(V![iv].2)

        x = limits[0][0] + x*griddeltas[0]
        y = limits[1][0] + y*griddeltas[1]
        z = limits[2][0] + z*griddeltas[2]

        VERTICES.append(Vector([x, y, z]))
        NORMALS.append(Vector([Double(N![iv].0), Double(N![iv].1), Double(N![iv].2)]))
    }

 
    var FACES = [[Int]]()

    for ie in 0..<Int(nT) {
        FACES.append([Int(T![ie].0), Int(T![ie].1), Int(T![ie].2) ])
    }

   

    //print("returning \(VERTICES.count) vertices, \(NORMALS.count) normals, \(FACES.count) faces")
    
    return (VERTICES,NORMALS,FACES)

}

public func generateTriangulation( probes:[Probe], probeRadius:Double, gridspacing:Double, 
        densityDelta:Double, densityEpsilon:Double, isoLevel:Double, numthreads:Int, axis:AXES=AXES.Z, mingridchunk:Int=100) 
            throws -> ([Vector], [Vector], [[Int]]) {
           // -> throws ([[Double]], [[Int]])  {


    let pcoords = probes .map { $0.center }

    var linearcoords = [Double]()

    for coord in pcoords {
        linearcoords += coord.coords
    }


    var spans = [[Double]]()

    for j in 0..<3 {
        let axiscoords = (0..<probes.count) .map { linearcoords[3*$0 + j] }
        spans.append([axiscoords.min()! - 2.0*probeRadius, axiscoords.max()! + 2.0*probeRadius])
    }


    // work with grid vertices, one more than number of grid divisions
    let gridvertices = (0..<3) .map { Int(round((spans[$0][1] - spans[$0][0])/gridspacing)) + 1  }

    let griddeltas = (0..<3) .map { (spans[$0][1] - spans[$0][0])/Double(gridvertices[$0] - 1) }

    print("\ngridvertices : \(gridvertices)  griddeltas : \(griddeltas)")



    // allot probes to nthreads
    //
    

    var probesForThread = [[Probe]]() 

    let probeChunk = probes.count / numthreads

    for ithread in 0..<numthreads {
        if ithread < numthreads - 1 {
            probesForThread.append( Array(probes[(ithread*probeChunk)..<((ithread+1)*probeChunk)]))
        }
        else {
            probesForThread.append( Array(probes[(ithread*probeChunk)..<probes.count]))
        }
    }

    // Get densities for probes by thread 

    // introduce major change - density blocks now will use X-Y-Z coordinate ordering, so do not need to do 
    // time-consuming operation of mapping every probe's local grid coordinates 
    // 
    // add into grid with normal expectation (Z-coord changing fastest), then use new transpose operation 
    // introduced into MathTools

    

    let group = DispatchGroup() 

    let gridShape = [ gridvertices[0], gridvertices[1], gridvertices[2] ]

    let gridDensity = Matrix<Double>(gridShape)

    let gridStrides = gridDensity.getStrides()


    for tidx in 0..<numthreads {

                print("enter thread \(tidx)")
       
                group.enter()

                computeQueue.async {

                        var data:[([Double],[Int])]?

                        let time0 = Date().timeIntervalSince1970

                        
                        let density =  densityForProbes( probes:probesForThread[tidx], radius:probeRadius, delta:densityDelta, 
                            epsilon:densityEpsilon, griddeltas:griddeltas, limits:spans, gridShape:gridShape, gridStrides:gridStrides )
                        
                        
                        blocksQueue.sync {
                            print("merge density vector from thread \(tidx)")

                            (0..<gridDensity.storage.count) .map { gridDensity.storage[$0] += density[$0]}
                        }

                        let time1 = Date().timeIntervalSince1970
                        print("leave thread \(tidx), time = \(time1 - time0), number probes = \(probesForThread[tidx].count)")
                        // adding this for some debugging info - something suddenly went wrong for no apparent reason
                        /*
                        let tots = data! .map { $0.0 . reduce( 0.0, + ) }
                        let totdens = tots .reduce( 0.0, + )
                        print("total density accumulated by thread \(tidx) = \(totdens)")

                        if totdens == 0.0 {
                            print("WARNING, total density for thread \(tidx) is zero ")
                            print("Probes for this thread : ")
                            for probe in probesForThread[tidx] {
                                var str = ""
                                str += "\(probe.center.coords[0]) , \(probe.center.coords[1]) , \(probe.center.coords[2])"
                                for atom in probe.atoms {
                                    str += " \(atom)"
                                }
                                print( str )
                            }
                        }
                        */
                        
                        group.leave()

                        
                    }

            
    } 

    group.wait()

    // Add probe densities into grid 

    // NOTE, for sanity's sake, I will feed marching cubes the transpose, so I have Z, Y, X densities;
    // 
    // Other arguments should be OK

    let time0 = Date().timeIntervalSince1970

    let gridDensity_T = gridDensity.transpose()

    let time1 = Date().timeIntervalSince1970

    print("\ntime for density matrix transpose = \(time1 - time0)")

    print("\nenter marching cubes ...")

    // As experiment, do not multithread, don't think I gain much once I need to stitch chunks together

    let tri = runMarchingCubes( density:gridDensity_T, limits:spans, 
                        griddeltas:griddeltas, gridvertices:gridvertices, isoLevel:isoLevel )



    return tri

}

 