/*
Copyright 2016-2018 Erik Erlandson

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

package org.isarnproject.sketches.java;

import java.lang.System;
import java.lang.StringBuilder;
import java.util.Arrays;
import java.util.Comparator;
import java.io.Serializable;
import java.util.concurrent.ThreadLocalRandom;

public final class TDigest implements Serializable {
    protected final double C;
    protected final int maxDiscrete;
    protected int nclusters = 0;
    protected double M = 0.0;
    protected double[] cent = null;
    protected double[] mass = null;
    protected double[] ftre = null;

    public TDigest() {
        this(COMPRESSION_DEFAULT, 0, INIT_SIZE_DEFAULT);
    }

    public TDigest(double compression) {
        this(compression, 0, INIT_SIZE_DEFAULT);
    }

    public TDigest(double compression, int maxDiscrete) {
        this(compression, maxDiscrete, INIT_SIZE_DEFAULT);
    }

    public TDigest(double compression, int maxDiscrete, int sz) {
        assert compression > 0.0;
        assert maxDiscrete >= 0;
        assert sz > 0;
        C = compression;
        this.maxDiscrete = maxDiscrete;
        cent = new double[sz];
        mass = new double[sz];
        ftre = new double[1 + sz];
        // ftre is 1-based. set ftre[0] to zero just to be tidy
        ftre[0] = 0.0;        
    }

    public TDigest(TDigest that) {
        C = that.C;
        maxDiscrete = that.maxDiscrete;
        nclusters = that.nclusters;
        M = that.M;
        cent = Arrays.copyOf(that.cent, nclusters);
        mass = Arrays.copyOf(that.mass, nclusters);
        ftre = Arrays.copyOf(that.ftre, nclusters);
    }

    public final void update(double x) {
        update(x, 1.0);
    }

    public final void update(double x, double w) {
        updateLogic(x, w);
        if ((nclusters > maxDiscrete) && (nclusters > R())) recluster();
    }

    private final void updateLogic(double x, double w) {
        if (nclusters == 0) {
            // clusters are empty, so (x,w) becomes the first cluster
            cent[0] = x;
            M = w;
            mass[0] = w;
            ftre[1] = w;
            nclusters += 1;
            return;
        }
        if (nclusters <= maxDiscrete) {
            // we are under the limit for discrete values to track
            int j = Arrays.binarySearch(cent, 0, nclusters, x);
            if (j >= 0) {
                // landed on existing cluster: add its mass and we're done
                M += w;
                mass[j] += w;
                ftInc(j, w);
            } else {
                // a new x value: insert as a new discrete cluster
                newCluster(-(j + 1), x, w);
            }
            return;
        }
        // get the index of the cluster closest to x
        int j = closest(x);
        if (x == cent[j]) {
            // landed on existing cluster: add its mass and we're done
            M += w;
            mass[j] += w;
            ftInc(j, w);
            return;
        }
        double m = mass[j];
        // q is the quantile of the closest cluster to x
        // (ftSum does the right thing (return 0) for j = 0)
        double q = (ftSum(j - 1) + (m / 2.0)) / M;
        // this is the upper-bound for the mass of closest cluster
        double ub = C * M * q * (1.0 - q);
        // dm is how much mass we're allowed to add to closest cluster
        double dm = Math.min(w, Math.max(0.0, ub - m));
        // rm is the remainder of the mass
        double rm = w - dm;
        if (dm > 0.0) {
            // Add any allowable mass to closest cluster and update its center.
            // It is safe to update center this way because it will remain
            // between x and original center, and so cannot move out of its original
            // ordering relative to its neighbors, because x is by previous logic
            // closer to cent[j] than any other cluster.
            double dc = dm * (x - cent[j]) / (m + dm);
            cent[j] += dc;
            M += dm;
            mass[j] += dm;
            ftInc(j, dm);
        }
        // if there is remaining mass, it becomes a new cluster
        if (rm > 0.0) newCluster((x < cent[j]) ? j : j + 1, x, rm);
    }

    public final void merge(TDigest that) {
        Integer[] indexes = new Integer[that.nclusters];
        for (int j = 0; j < that.nclusters; ++j) indexes[j] = j;
        // sort so that largest clusters are first.
        // inserting large to small yields stable distribution estimations
        Comparator<Integer> cmp = new Comparator<Integer>() {
            @Override
            public int compare(Integer a, Integer b) {
                return (int)Math.signum(that.mass[b] - that.mass[a]);
            }
        };
        Arrays.sort(indexes, cmp);
        for (int j: indexes) update(that.cent[j], that.mass[j]);
    }

    public final void recluster() {
        int[] indexes = new int[nclusters];
        for (int j = 0; j < nclusters; ++j) indexes[j] = j;
        intShuffle(indexes);
        int sz = cent.length;
        double[] oldCent = cent;
        double[] oldMass = mass;
        cent = new double[sz];
        mass = new double[sz];
        reset();
        for (int j: indexes) updateLogic(oldCent[j], oldMass[j]);
    }

    public final void reset() {
        nclusters = 0;
        M = 0.0;
    }

    private final void newCluster(int j, double x, double w) {
        double[] newCent = cent;
        double[] newMass = mass;
        double[] newFtre = ftre;
        int sz = cent.length;
        if (nclusters >= sz) {
            int szinc = (int)Math.ceil(0.1 * (double)sz);
            sz += szinc;
            newCent = new double[sz];
            newMass = new double[sz];
            newFtre = new double[1 + sz];
            System.arraycopy(cent, 0, newCent, 0, j);
            System.arraycopy(mass, 0, newMass, 0, j);
        }
        // arraycopy can handle when cent == newCent
        System.arraycopy(cent, j, newCent, 1 + j, nclusters - j);
        System.arraycopy(mass, j, newMass, 1 + j, nclusters - j);
        // do this after copies above
        newCent[j] = x;
        newMass[j] = w;
        nclusters += 1;
        cent = newCent;
        mass = newMass;
        ftre = newFtre;
        Arrays.fill(ftre, 0, 1 + nclusters, 0.0);
        for (int k = 0; k < nclusters; ++k) ftInc(k, mass[k]);
        M += w;
    }

    private final int closest(double x) {
        int j = Arrays.binarySearch(cent, 0, nclusters, x);
        // exact match, return its index:
        if (j >= 0) return j;
        // x is not a cluster center, get its insertion index:
        j = -(j + 1);
        // x is to left of left-most cluster:
        if (j == 0) return j;
        // x is to right of right-most cluster:
        if (j == nclusters) return j - 1;
        // x is between two clusters, return index of closest:
        double dL = x - cent[j - 1];
        double dR = cent[j] - x;
        return (dL < dR) ? (j - 1) : j;
    }

    public final int size() {
        return nclusters;
    }

    public final double mass() {
        return M;
    }

    public final double getCompression() {
        return C;
    }

    public final int getMaxDiscrete() {
        return maxDiscrete;
    }

    public final double[] getCent() {
        return cent;
    }

    public final double[] getMass() {
        return mass;
    }

    public final double[] getFT() {
        return ftre;
    }
    
    public final boolean isEmpty() {
        return nclusters == 0;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder("TDigest(");
        for (int j = 0; j < nclusters; ++j) {
            if (j > 25) {
                sb.append(" ...");
                break;
            }
            if (j > 0) sb.append(", ");
            sb.append(cent[j])
                .append(" -> (")
                .append(mass[j])
                .append(", ")
                .append(ftSum(j))
                .append(")");
        }
        sb.append(")");
        return sb.toString();
    }

    public final double samplePDF() {
        return cdfInverse(ThreadLocalRandom.current().nextDouble());
    }

    public final double samplePMF() {
        return cdfDiscreteInverse(ThreadLocalRandom.current().nextDouble());
    }

    public final double sample() {
        if (nclusters <= maxDiscrete) {
            return cdfDiscreteInverse(ThreadLocalRandom.current().nextDouble());
        } else {
            return cdfInverse(ThreadLocalRandom.current().nextDouble());
        }    
    }

    public final double cdf(double x) {
        int j1 = rcovj(x);
        if (j1 < 0) return 0.0;
        if (j1 >= nclusters - 1) return 1.0;
        int j2 = j1 + 1;
        double c1 = cent[j1];
        double c2 = cent[j2];
        double tm1 = mass[j1];
        double tm2 = mass[j2];
        double s = ftSum(j1 - 1);
        double d1 = (j1 == 0) ? 0.0 : tm1 / 2.0;
        double m1 = s + d1;
        double m2 = m1 + (tm1 - d1) + ((j2 == nclusters - 1) ? tm2 : tm2 / 2.0);
        double m = m1 + (x - c1) * (m2 - m1) / (c2 - c1);
        return Math.min(m2, Math.max(m1, m)) / M;
    }

    public final double cdfDiscrete(double x) {
        int j = rcovj(x);
        return ftSum(j) / M;
    }

    public final double cdfInverse(double q) {
        if (q < 0.0 || q > 1.0) return Double.NaN;
        if (nclusters == 0) return Double.NaN;
        if (nclusters == 1) return cent[0];
        double m = q * M;
        int j1 = rmcovj(m);
        int j2 = j1 + 1;
        double c1 = cent[j1];
        double c2 = cent[j2];
        double tm1 = mass[j1];
        double tm2 = mass[j2];
        double s = ftSum(j1 - 1);
        double d1 = (j1 == 0) ? 0.0 : tm1 / 2.0;
        double m1 = s + d1;
        double m2 = m1 + (tm1 - d1) + ((j2 == nclusters - 1) ? tm2 : tm2 / 2.0);
        double x = c1 + (m - m1) * (c2 - c1) / (m2 - m1);
        return Math.min(c2, Math.max(c1, x));
    }

    public final double cdfDiscreteInverse(double q) {
        if (q < 0.0 || q > 1.0) return Double.NaN;
        if (nclusters == 0) return Double.NaN;
        if (nclusters == 1) return cent[0];
        double m = q * M;
        int j = lmcovj(m);
        return cent[j];
    }

    // returns the index of a right mass cover
    // ftSum(j) <= m < ftSum(j+1)
    private final int rmcovj(double m) {
        assert nclusters >= 2;
        assert (m >= 0.0) && (m <= M);
        int beg = 0;
        double mbeg = 0.0;
        int end = nclusters - 1;
        double mend = M;
        while ((end - beg) > 1) {
            int mid = (beg + end) / 2;
            double mmid = ftSum(mid);
            if (m >= mmid) {
                beg = mid;
                mbeg = mmid;
            } else {
                end = mid;
                mend = mmid;
            }
        }
        return beg;
    }

    // returns the index of a left mass cover
    // ftSum(j-1) < m <= ftSum(j)
    private final int lmcovj(double m) {
        assert nclusters >= 2;
        assert (m >= 0.0) && (m <= M);
        int beg = -1;
        double mbeg = 0.0;
        int end = nclusters - 1;
        double mend = M;
        while ((end - beg) > 1) {
            int mid = (beg + end) / 2;
            double mmid = ftSum(mid);
            if (m <= mmid) {
                end = mid;
                mend = mmid;
            } else {
                beg = mid;
                mbeg = mmid;
            }
        }
        return end;
    }

    // returns the left index of a right-cover
    private final int rcovj(double x) {
        int j = Arrays.binarySearch(cent, 0, nclusters, x);
        // exact match, return its index:
        if (j >= 0) return j;
        // x is not a cluster center, get its insertion index:
        j = -(j + 1);
        // x is to left of left-most cluster:
        if (j == 0) return -1;
        // return the index to the left of x:
        return j - 1;
    }

    // cumulative-sum algorithm for a Fenwick tree
    private final double ftSum(int j) {
        j += 1;
        double s = 0.0;
        while (j > 0) {
            s += ftre[j];
            j -= j & (-j); // dec by least significant nonzero bit of j
        }
        return s;
    }

    // increment algorithm for a Fenwick tree
    private final void ftInc(int j, double w) {
        j += 1;
        while (j <= nclusters) {
            ftre[j] += w;
            j += j & (-j); // inc by least significant nonzero bit of j
        }
    }

    @Override
    public boolean equals(Object that) {
        if (!(that instanceof TDigest)) return false;
        if (this == that) return true;
        TDigest rhs = (TDigest)that;
        if (C != rhs.C) return false;
        if (maxDiscrete != rhs.maxDiscrete) return false;
        if (nclusters != rhs.nclusters) return false;
        if (M != rhs.M) return false;
        if (!equal(cent, rhs.cent, nclusters)) return false;
        if (!equal(mass, rhs.mass, nclusters)) return false;
        // if masses are equal, cumulative ftre had better also be equal
        return true;
    }

    // I can't believe java just added this to Arrays in java 9
    static final boolean equal(double[] lhs, double[] rhs, int n) {
        for (int j = 0; j < n; ++j) {
            if (lhs[j] != rhs[j]) return false;
        }
        return true;
    }

    @Override
    public int hashCode() {
        int h = nclusters;
        h ^= doubleHash(M);
        if (nclusters >= 1) {
            h ^= doubleHash(cent[0]);
            h ^= doubleHash(mass[0]);
            h ^= doubleHash(ftre[1]);
        }
        if (nclusters >= 2) {
            h ^= doubleHash(cent[nclusters - 1]);
            h ^= doubleHash(mass[nclusters - 1]);
            h ^= doubleHash(ftre[nclusters]);
        }
        if (nclusters >= 3) {
            int j = nclusters / 2;
            h ^= doubleHash(cent[j]);
            h ^= doubleHash(mass[j]);
            h ^= doubleHash(ftre[1 + j]);
        }
        return h;
    }

    // I can't believe Double doesn't provide a static method for this
    static final int doubleHash(double x) {
        long v = Double.doubleToLongBits(x);
        return (int)(v ^ (v >>> 32));
    }

    protected final int R() {
        return (int)(K / C);
    }

    public static final double K = 10.0 * 50.0;

    /** delta * E[clusters] ~ 50 */
    public static final double COMPRESSION_DEFAULT = 50.0 / 100.0;

    public static final int INIT_SIZE_DEFAULT = 5;

    public static TDigest empty() {
        return new TDigest(COMPRESSION_DEFAULT, 0, INIT_SIZE_DEFAULT);
    }

    public static TDigest empty(double compression) {
        return new TDigest(compression, 0, INIT_SIZE_DEFAULT);
    }

    public static TDigest empty(double compression, int maxDiscrete) {
        return new TDigest(compression, maxDiscrete, INIT_SIZE_DEFAULT);
    }

    public static TDigest empty(double compression, int maxDiscrete, int sz) {
        return new TDigest(compression, maxDiscrete, sz);
    }

    public static TDigest merge(TDigest ltd, TDigest rtd) {
        if (ltd.size() < rtd.size()) return merge(rtd, ltd);
        if (rtd.size() == 0) return ltd;
        if (rtd.size() == 1) {
            ltd.update(rtd.cent[0], rtd.mass[0]);
            return ltd;
        }
        if (rtd.mass() < ltd.mass()) {
            ltd.merge(rtd);
            return ltd;
        } else {
            rtd.merge(ltd);
            return rtd;
        }
    }

    public static TDigest sketch(double[] data) {
        return sketch(data, COMPRESSION_DEFAULT, 0, INIT_SIZE_DEFAULT);
    }

    public static TDigest sketch(double[] data, double compression) {
        return sketch(data, compression, 0, INIT_SIZE_DEFAULT);
    }

    public static TDigest sketch(double[] data, double compression, int maxDiscrete) {
        return sketch(data, compression, maxDiscrete, INIT_SIZE_DEFAULT);
    }

    public static TDigest sketch(double[] data, double compression, int maxDiscrete, int sz) {
        TDigest td = empty(compression, maxDiscrete, sz);
        for (double x: data) td.update(x, 1.0);
        if (td.size() > maxDiscrete) td.recluster();
        return td;
    }

    public static void intShuffle(int[] data) {
        intShuffle(data, 0, data.length);
    }

    public static void intShuffle(int[] data, int end) {
        intShuffle(data, 0, end);
    }

    public static void intShuffle(int[] data, int beg, int end) {
        ThreadLocalRandom rnd = ThreadLocalRandom.current();
        end -= 1;
        while (end > beg) {
            int r = rnd.nextInt(beg, end);
            int d = data[end];
            data[end] = data[r];
            data[r] = d;
            end -= 1;
        }
    }
}
