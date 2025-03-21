#pragma once

#include <stdio.h>
#include <string>
#include <vector>
#include <filesystem>
#include <stack>
#include <memory>
#include "mybasics.h"

struct Indent
{
    Indent(struct ISerializer& s, const char* sName) : m_serializer(s) { }
    ~Indent();
private:
    ISerializer& m_serializer;
};

struct ISerializer
{
    virtual ~ISerializer()
    {

    }

    template <class T>
    void serializeArraySize(const char *sName, std::vector<T>& p)
    {
        size_t size = p.size();
        serializeSimpleType("arraySize", size);
        p.resize(size);
    }

    template <class T>
    void serializeStdArray(const char *sName, std::vector<T>& p)
    {
        std::shared_ptr<Indent> pIndent = pushIndent(sName);
        serializeArraySize(sName, p);
        if (p.size() > 0)
        {
            serializePreallocatedMem(sName, &p[0], sizeof(p[0]) * (NvU32)p.size());
        }
    }

    template <class T>
    void serializeArrayOfPointers(const char *sName, std::vector<T*>& p)
    {
        std::shared_ptr<Indent> pIndent = pushIndent(sName);
        // serialize array size
        serializeArraySize(sName, p);

        // serialize each pointer in the array
        for (NvU32 u = 0; u < p.size(); ++u)
        {
            NvU32 ptrExists = (p[u] != nullptr);
            serializeSimpleType("ptrExists", ptrExists);
            if (!ptrExists)
            {
                continue;
            }
            if (m_isReading)
            {
                p[u] = new T();
            }
            char sBuffer[16];
            sprintf_s(sBuffer, "[%d]", u);
            std::shared_ptr<Indent> pIndent = pushIndent(sBuffer);
            p[u]->serialize(sName, *this);
        }
    }

    template <class T>
    void serializeSharedPtr(const char *sName, std::shared_ptr<T>& ref)
    {
        NvU32 refExists = (ref != nullptr);
        serializeSimpleType("refExists", refExists);
        if (!refExists)
        {
            return;
        }
        if (ref == nullptr)
        {
            ref = std::make_shared<T>();
        }
        ref->serialize(sName, *this);
    }

    template <class T>
    void serializeArrayOfSharedPtrs(const char *sName, std::vector<std::shared_ptr<T>>& refs)
    {
        std::shared_ptr<Indent> pIndent = pushIndent(sName);
        serializeArraySize(sName, refs);
        
        // serialize each shared_ptr in the array
        for (NvU32 u = 0; u < refs.size(); ++u)
        {
            char sBuffer[16];
            sprintf_s(sBuffer, "[%d]", u);
            std::shared_ptr<Indent> pIndent = pushIndent(sBuffer);
            serializeSharedPtr(nullptr , refs[u]);
        }
    }

    template <class T>
    void serializeSimpleType(const char *sName, T & value)
    {
        serializePreallocatedMem(sName, &value, sizeof(value));
    }
    virtual void serializePreallocatedMem(const char *sName, void* pMem, NvU32 memSizeInBytes) = 0;

    virtual std::shared_ptr<Indent> pushIndent(const char* sName)
    {
        return nullptr; // by default indents are ignored
    }
    virtual void notifyIndentDestroyed(Indent *p) { }

protected:
    bool m_isReading = false;
    std::wstring makeFullPath(const std::filesystem::path& path)
    {
        if (path.is_absolute())
            return path;
        return std::wstring(L"c:\\atomNets\\") + path.c_str();
    }
};

struct MyWriter : public ISerializer
{
    MyWriter(const std::filesystem::path& path)
    {
        std::wstring sFullPath = makeFullPath(path);
        _wfopen_s(&m_fp, sFullPath.c_str(), L"wb");
        nvAssert(m_fp != nullptr);
    }
    ~MyWriter()
    {
        fclose(m_fp);
    }
    virtual void serializePreallocatedMem(const char *sName, void* pMem, NvU32 memSizeInBytes) override
    {
        if (memSizeInBytes > 0)
        {
            nvAssert(m_fp != nullptr);
            fwrite(pMem, 1, memSizeInBytes, m_fp);
        }
    }
private:
    FILE* m_fp = nullptr;
};

struct MyReader : public ISerializer
{
    MyReader(const std::filesystem::path& path)
    {
        std::wstring sFullPath = makeFullPath(path);
        _wfopen_s(&m_fp, sFullPath.c_str(), L"rb");
        nvAssert(m_fp != nullptr);
        m_isReading = true;
    }
    ~MyReader()
    {
        fclose(m_fp);
    }
    virtual void serializePreallocatedMem(const char *sName, void* pMem, NvU32 memSizeInBytes) override
    {
        nvAssert(m_fp != nullptr);
        fread(pMem, 1, memSizeInBytes, m_fp);
    }
private:
    FILE* m_fp = nullptr;
};

struct TextWriter : public ISerializer
{
    TextWriter(char* sFileName)
    {
        fopen_s(&m_fp, sFileName, "wt");
    }
    ~TextWriter()
    {
        nvAssert(m_fp != nullptr);
        fclose(m_fp);
    }

    virtual void serializePreallocatedMem(const char *sName, void* pMem, NvU32 memSizeInBytes) override
    {
        printCurIndent();
        nvAssert(m_fp != nullptr);
        fprintf(m_fp, "%s:\n", sName);
        if (memSizeInBytes % sizeof(float) == 0)
        {
            int n = memSizeInBytes / sizeof(float);
            printCurIndent();
            fprintf(m_fp, "  ");
            float* pFloats = (float*)pMem;
            if (hasDenorms(pFloats, n))
            {
                int* pInts = (int*)pMem;
                for (int i = 0; i < n; ++i)
                {
                    fprintf(m_fp, "%d, ", pInts[i]);
                }
            }
            else
            {
                for (int i = 0; i < n; ++i)
                {
                    fprintf(m_fp, "%f, ", pFloats[i]);
                }
            }
            fprintf(m_fp, "\n");
        }
    }

    virtual std::shared_ptr<Indent> pushIndent(const char* sName) override
    {
        printCurIndent();
        nvAssert(m_fp != nullptr);
        fprintf(m_fp, "%s {\n", sName);
        std::shared_ptr<Indent> p = std::make_shared<Indent>(*this, sName);
        m_pIndents.push(p.get());
        return p;
    }
    virtual void notifyIndentDestroyed(Indent* p) override
    {
        nvAssert(p == m_pIndents.top());
        m_pIndents.pop();
        printCurIndent();
        fprintf(m_fp, "}\n");
    }

private:
    bool hasDenorms(float* pFloats, NvU32 n)
    {
        int* pInts = (int*)pFloats;
        n = std::min(n, 30u);
        for (NvU32 u = 0; u < n; ++u)
        {
            if (pInts[u] > 0 && pInts[u] < 30)
                return true;
        }
        return false;
    }
    void printCurIndent()
    {
        nvAssert(m_fp != nullptr);
        for (int i = 0; i < m_pIndents.size(); ++i)
        {
            fprintf(m_fp, "  ");
        }
    }
    FILE* m_fp = nullptr;
    std::stack<Indent*> m_pIndents;
};

inline Indent::~Indent()
{
    m_serializer.notifyIndentDestroyed(this);
}